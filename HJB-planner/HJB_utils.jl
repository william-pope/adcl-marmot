# HJB_utils.jl

using StaticArrays

include("dynamics_models.jl")

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
    E_std = [[-veh.ext2axle, 1/2*veh.ext_w],         # top left
            [-veh.ext2axle, -1/2*veh.ext_w],         # bottom left
            [-veh.ext2axle+veh.ext_l, -1/2*veh.ext_w],   # bottom right
            [-veh.ext2axle+veh.ext_l, 1/2*veh.ext_w],    # top right
            [-veh.ext2axle, 1/6*veh.ext_w],              # left edge refinement              
            [-veh.ext2axle, -1/6*veh.ext_w],
            [-veh.ext2axle+1/4*veh.ext_l, -1/2*veh.ext_w],   # bottom edge refinement
            [-veh.ext2axle+1/2*veh.ext_l, -1/2*veh.ext_w],
            [-veh.ext2axle+3/4*veh.ext_l, -1/2*veh.ext_w],
            [-veh.ext2axle+veh.ext_l, 1/6*veh.ext_w],        # right edge refinement
            [-veh.ext2axle+veh.ext_l, -1/6*veh.ext_w],
            [-veh.ext2axle+1/4*veh.ext_l, 1/2*veh.ext_w],    # top edge refinement
            [-veh.ext2axle+1/2*veh.ext_l, 1/2*veh.ext_w],
            [-veh.ext2axle+3/4*veh.ext_l, 1/2*veh.ext_w]]

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

function find_idx(val, array)
    if val in array
        idx = indexin(val, array)[1]    # SPEED: allocation (think it creates an array sometimes -> bad)
                                        #   - 101.8 μs when array, 8.7 μs when no array (in practice will never exactly be on grid node)
    else
        idx = searchsortedfirst(array, val) - 1
    end

    return idx
end

function plot_polygon(my_plot, P, lw, lc, ll)
    P_x_pts = [P[:,1]; P[1,1]]
    P_y_pts = [P[:,2]; P[1,2]]

    plot!(my_plot, P_x_pts, P_y_pts, linewidth=lw, linecolor=lc, label=ll)
end

# # TO-DO: make sure this isn't causing any issues, conflicts with other boundary mechanisms
#     # adjust position when near edge for boundary issues
#     if abs(x_k[1] - env.x_grid[1]) < 2*env.h_xy
#         x_k[1] = env.x_grid[3]
#     elseif abs(x_k[1] - env.x_grid[end]) < 2*env.h_xy
#         x_k[1] = env.x_grid[end-2]
#     end

#     if abs(x_k[2] - env.y_grid[1]) < 2*env.h_xy
#         x_k[2] = env.y_grid[3]
#     elseif abs(x_k[2] - env.y_grid[end]) < 2*env.h_xy
#         x_k[2] = env.y_grid[end-2]
#     end

#     x_k[3] > pi ? x_k[3] = x_k[3] - 2*pi : x_k[3] = x_k[3]

#     println("\n--- ---\nx_k: ", x_k)

# # ISSUE: dealing with obstacles during gradient approximation
# #   - unsure what to do when one of interpolation nodes is in an obstacle
# #   - depending on grid size, interp node could be somewhat far away from actual point of interest
# #   - need to somehow calculate a gradient without using affected nodes
# #   - gradient process is only a function of the state, not related to actions being tested
# #   - obstacles are not a part of the flow field, treat like blank regions
# #   - partials can be set to 0 if needed, just implies a flat gradient along that axis
# #       - probably more accurate to switch to single point difference method, but will be a little harder

# function value_gradient(x::Vector{Float64}, U::Array{Float64, 3}, O::Array{Bool, 3}, env::Environment)
#     # println(x)

#     # ISSUE: chattering still happening near obstacles
#     #   - in this case, think neighbor in x and theta are obstacle states
#     #   - actually negative y is too (doesn't look like it should...)
#     #   - obstacles on all three axes -> frozen
#     #   - single point gradient will fix this

#     #   - each axis gradient is independent, can look at each in a for loop
#     #   - logic:
#     #       - if both free, use center difference method with both points
#     #       - if one point in obs, use single diff method with other point and x
#     #       - if both in obs, set du_dy = 0.0

#     #   - think issue is due to both outlier point and central point touching obstacle
    
#     # interpolate value at surrounding points
#     v_x, obs_x = value_interp(x, U, O, env)
#     v_1p, obs_1p = value_interp(x + [env.h_xy, 0.0, 0.0], U, O, env)
#     v_1n, obs_1n = value_interp(x - [env.h_xy, 0.0, 0.0], U, O, env)
#     v_2p, obs_2p = value_interp(x + [0.0, env.h_xy, 0.0], U, O, env)
#     v_2n, obs_2n = value_interp(x - [0.0, env.h_xy, 0.0], U, O, env)
#     v_3p, obs_3p = value_interp(x + [0.0, 0.0, env.h_theta], U, O, env)
#     v_3n, obs_3n = value_interp(x - [0.0, 0.0, env.h_theta], U, O, env)

#     # ISSUE: interpolation at center point can be in obstacle
#     #   - occurs just after vehicle rounds corner of circle (9 o'clock position)
#     #   - y_neg in obstacle seems strange, since behind vehicle
#     #   - this might just be what happens on a diagonal face
#     #   - happens when making sharp turns around obstacles (need more theta steps?)
#     #   - is step size just too large? seems like we're already at the limit

#     # calculating du_dx
#     if obs_1p == false && obs_1n == false
#         du_dx = (v_1p - v_1n)/(2*env.h_xy)
#     elseif obs_1p == true && obs_1n == false
#         du_dx = (v_x - v_1n)/env.h_xy
#     elseif obs_1p == false && obs_1n == true
#         du_dx = (v_1p - v_x)/env.h_xy
#     else
#         du_dx = 0.0
#     end

#     # calculating du_dy
#     if obs_2p == false && obs_2n == false
#         du_dy = (v_2p - v_2n)/(2*env.h_xy)
#     elseif obs_2p == true && obs_2n == false
#         du_dy = (v_x - v_2n)/env.h_xy
#     elseif obs_2p == false && obs_2n == true
#         du_dy = (v_2p - v_x)/env.h_xy
#     else
#         du_dy = 0.0
#     end

#     # calculating du_dtheta
#     if obs_3p == false && obs_3n == false
#         du_dtheta = (v_3p - v_3n)/(2*env.h_theta)
#     elseif obs_3p == true && obs_3n == false
#         du_dtheta = (v_x - v_3n)/env.h_theta
#     elseif obs_3p == false && obs_3n == true
#         du_dtheta = (v_3p - v_x)/env.h_theta
#     else
#         du_dtheta = 0.0
#     end

#     dV = [du_dx, du_dy, du_dtheta]

#     # println("dV: ", dV)

#     return dV
# end

#     # # calculate gradient at current state
#     # dV = value_gradient(x_k, U, O, env)       # SPEED: big focus

#     # vdot_min = Inf
#     # a_opt = A[1]
#     # for i in 1:size(A,1)
#     #     xs = SVector{3, Float64}(x_k)     # TO-DO: clean up variable types
#     #     xdot = car_EoM(xs, A[i], veh)

#     #     vdot = dot(dV, xdot)

#     #     if vdot < vdot_min
#     #         vdot_min = vdot
#     #         a_opt = A[i]
#     #     end
#     # end