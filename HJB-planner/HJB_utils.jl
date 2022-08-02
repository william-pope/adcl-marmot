# HJB_utils.jl

using StaticArrays

# define discretized grid struct
struct Environment
    h_xy::Float64
    h_theta::Float64
    x_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    y_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    theta_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    W::Matrix{Float64}
    T_xy::Matrix{Float64}
    T_theta::Vector{Vector{Float64}}
    O_vec::Vector{Matrix{Float64}}
end

struct Vehicle
    c_vf::Float64   # max forward speed [m/s]
    c_vb::Float64   # max backward speed [m/s]
    c_phi::Float64  # max steering angle [rad]
    wb::Float64     # wheelbase [m]
    l::Float64      # length [m]
    w::Float64      # width [m]
    b2a::Float64    # rear bumber to rear axle [m]
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k::Vector{Float64}, u::Vector{Float64}, dt, EoM::Function, veh::Vehicle)
    x_ks = SVector{3, Float64}(x_k)

    w1 = EoM(x_ks, u, veh)
    w2 = EoM(x_ks + w1*dt/2, u, veh)
    w3 = EoM(x_ks + w2*dt/2, u, veh)
    w4 = EoM(x_ks + w3*dt, u, veh)

    x_k1s = x_ks + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    # convert x_k1 back
    x_k1 = [x_k1s[1], x_k1s[2], x_k1s[3]]   # TO-DO: clean up variable types

    # adjust theta within bounds
    if x_k1[3] > pi
        x_k1[3] -= 2*pi
    elseif x_k1[3] < -pi
        x_k1[3] += 2*pi
    end

    return [x_k1[1], x_k1[2], x_k1[3]]
end

# target set checker
function in_target_set(y, env::Environment, veh::Vehicle)
    # check position of corners
    C = pose_to_edges(y, veh)
    rows = [collect(1:size(env.T_xy, 1)); 1]    # SPEED: allocation

    for c in C
        for i in 1:size(env.T_xy, 1)
            i1 = rows[i]
            i2 = rows[i+1]

            x1 = env.T_xy[i1,1]
            y1 = env.T_xy[i1,2]
            x2 = env.T_xy[i2,1]
            y2 = env.T_xy[i2,2]

            val = (y1 - y2)*c[1] + (x2 - x1)*c[2] + x1*y2 - x2*y1

            if val < 0
                return false
            end
        end
    end

    # check orientation
    in_theta = zeros(Bool, size(env.T_theta, 1))
    for i in 1:size(env.T_theta, 1)
        if y[3] >= minimum(env.T_theta[i]) && y[3] <= maximum(env.T_theta[i])
            in_theta[i] = 1
        end
    end
    if any(in_theta) == false
        return false
    end

    return true
end

# obstacle set checker
function in_obstacle_set(y, env::Environment, veh::Vehicle)
    E = pose_to_edges(y, veh::Vehicle)

    # add more points "c" along sides of vehicle
    for e in E 
        for Oi in env.O_vec
            rows = [collect(1:size(Oi, 1)); 1]  # TO-DO: fairly large impact on runtime

            ineqs = zeros(Bool, size(Oi, 1))
            for i in 1:size(Oi, 1)
                i1 = rows[i]
                i2 = rows[i+1]

                x1 = Oi[i1,1]
                y1 = Oi[i1,2]
                x2 = Oi[i2,1]
                y2 = Oi[i2,2]

                val = (y1 - y2)*e[1] + (x2 - x1)*e[2] + x1*y2 - x2*y1

                if val >= 0
                    ineqs[i] = 1
                end
            end

            if all(ineqs) == true
                return true
            end
        end
    end

    return false
end

function circle_to_polygon(OC_cir)
    x_c = OC_cir[1]
    y_c = OC_cir[2]
    r_c = OC_cir[3]

    # number of points used to discretize edge of circle
    pts = 8

    # circle radius is used as midpoint radius for polygon faces
    r_p = r_c/cos(pi/pts)

    theta_rng = range(0, 2*pi, pts+1)
    OC = Array{Float64}(undef, pts+1, 2)

    for (i, theta) in enumerate(theta_rng)
        OC[i,1] = x_c + r_p*cos(theta)
        OC[i,2] = y_c + r_p*sin(theta)
    end

    OC = OC[1:end-1, 1:end]

    return OC
end

# workspace checker
function in_workspace(y, env::Environment, veh::Vehicle)
    E = pose_to_edges(y, veh::Vehicle)

    for e in E 
        rows = [collect(1:size(env.W, 1)); 1]

        for i in 1:size(env.W, 1)
            i1 = rows[i]
            i2 = rows[i+1]

            x1 = env.W[i1,1]
            y1 = env.W[i1,2]
            x2 = env.W[i2,1]
            y2 = env.W[i2,2]

            val = (y1 - y2)*e[1] + (x2 - x1)*e[2] + x1*y2 - x2*y1

            if val <= 0
                return false
            end
        end
    end

    return true
end

function on_boundary(i, j, k, env::Environment)
    if i == 1 || i == size(env.x_grid,1) || j == 1 || j == size(env.y_grid,1)
        return true
    else
        return false
    end
end

# calculates position of corners of rectangular vehicle
function pose_to_edges(x, veh::Vehicle)
    # define standard points
    # SPEED: allocation, should be able to define as StaticArray
    E_std = [[-veh.b2a, 1/2*veh.w],         # top left
            [-veh.b2a, -1/2*veh.w],         # bottom left
            [-veh.b2a+veh.l, -1/2*veh.w],   # bottom right
            [-veh.b2a+veh.l, 1/2*veh.w],    # top right
            [-veh.b2a, 1/6*veh.w],              # left edge refinement              
            [-veh.b2a, -1/6*veh.w],
            [-veh.b2a+1/4*veh.l, -1/2*veh.w],   # bottom edge refinement
            [-veh.b2a+1/2*veh.l, -1/2*veh.w],
            [-veh.b2a+3/4*veh.l, -1/2*veh.w],
            [-veh.b2a+veh.l, 1/6*veh.w],        # right edge refinement
            [-veh.b2a+veh.l, -1/6*veh.w],
            [-veh.b2a+1/4*veh.l, 1/2*veh.w],    # top edge refinement
            [-veh.b2a+1/2*veh.l, 1/2*veh.w],
            [-veh.b2a+3/4*veh.l, 1/2*veh.w]]

    # calculate points for current pose
    E_rt = [Vector{Float64}(undef, 2) for _ = 1:size(E_std,1)]  # SPEED: allocation, might be fixed by SA above
    for (i, e_std) in enumerate(E_std)
        # E_rt[i] = RotMatrix{2}(x[3])*e_std + [x[1], x[2]]   # SPEED: allocation

        E_rt[i] = RotMatrix{2}(x[3])*e_std
        E_rt[i][1] += x[1]
        E_rt[i][2] += x[2]
    end

    return E_rt
end

function plot_polygon(my_plot, P, lw, lc, ll)
    P_x_pts = [P[:,1]; P[1,1]]
    P_y_pts = [P[:,2]; P[1,2]]

    plot!(my_plot, P_x_pts, P_y_pts, linewidth=lw, linecolor=lc, label=ll)
end