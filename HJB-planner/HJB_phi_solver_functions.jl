

using StaticArrays
using LinearAlgebra
using Rotations
using BenchmarkTools
using ProfileView

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"

# define discretized grid struct
struct Environment
    h_xy::Float64
    h_theta::Float64
    h_phi::Float64
    x_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    y_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    theta_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    phi_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    W::Matrix{Float64}
    T_xy::Matrix{Float64}
    T_theta::Vector{Vector{Float64}}
    O_vec::Vector{Matrix{Float64}}
end

struct Vehicle
    c_vf::Float64   # max forward speed [m/s]
    c_vb::Float64   # max backward speed [m/s]
    c_phi::Float64  # max steering angle [rad]
    c_phi_dot::Float64  # max steering angular velocity [rad/s]
    wb::Float64     # wheelbase [m]
    l::Float64      # length [m]
    w::Float64      # width [m]
    b2a::Float64    # rear bumber to rear axle [m]
end

# target set checker
function in_target_set(x, env::Environment, veh::Vehicle)
    # check position of corners
    C = pose_to_edges(x, veh)
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
        if x[3] >= minimum(env.T_theta[i]) && x[3] <= maximum(env.T_theta[i])
            in_theta[i] = 1
        end
    end
    if any(in_theta) == false
        return false
    end

    return true
end

# obstacle set checker
function in_obstacle_set(x, env::Environment, veh::Vehicle)
    E = pose_to_edges(x, veh::Vehicle)

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

# TO-DO: should enforce phi limits here?
#   - could rename to "in_state_space"
#   - look at where function is used, context may tell

# workspace checker
function in_workspace(x, env::Environment, veh::Vehicle)
    E = pose_to_edges(x, veh::Vehicle)

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

# (?): same purpose as in_workspace?
#   - used once, in initialize_array()
#   - used for same purpose as in_obstacle(), assigns bad value to grid node
function on_boundary(i, j, k, p, env::Environment)
    if (i == 1 || i == size(env.x_grid,1)) || (j == 1 || j == size(env.y_grid,1)) || (p == 1 || p == size(env.phi_grid,1))
        return true
    else
        return false
    end
end

# NOTE: want edges of phi grid to be included in state space, so that vehicle can use max steering
#   - ways to implement:
#       1) add extra node on each end of phi_grid, label as boundary (would need to enforce constraint on any point between the two)
#       2) treat phi_max as inclusive boundary, solve for nodes in HJB_solver (need to account for missing neighbor to outside)
#   - would need to enforce hard constraint past phi_max in either case
#       - (?): where does this get enforced? probably only in planning
#   - method 1) is closest to existing method

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

# initialize value approximations
function initialize_value_array(env::Environment, veh::Vehicle)
    Up = zeros(Float64, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.phi_grid,1)))
    Ip = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.phi_grid,1)))
    T = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.phi_grid,1)))
    O = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.phi_grid,1)))

    for i in 1:size(env.x_grid,1)
        for j in 1:size(env.y_grid,1)
            for k in 1:size(env.theta_grid,1)
                for p in 1:size(env.phi_grid,1)
                    x_i = env.x_grid[i]
                    y_j = env.y_grid[j]
                    theta_k = env.theta_grid[k]
                    phi_k = env.phi_grid[p]

                    x_ijkp = [x_i, y_j, theta_k, phi_k]

                    if in_target_set(x_ijkp, env, veh) == true
                        Up[i,j,k,p] = 0.0
                        Ip[i,j,k,p] = true
                        T[i,j,k,p] = true
                    elseif on_boundary(i, j, k, p, env) == true || in_obstacle_set(x_ijkp, env, veh) == true
                        Up[i,j,k,p] = 1000.0
                        Ip[i,j,k,p] = false
                        O[i,j,k,p] = true
                    else
                        Up[i,j,k,p] = 100.0
                        Ip[i,j,k,p] = false
                        T[i,j,k,p] = false
                        O[i,j,k,p] = false
                    end
                end
            end
        end
    end

    U = deepcopy(Up)
    I = deepcopy(Ip)

    return U, Up, I, Ip, T, O
end

# use modules
# Julia workflows

# @code_warntype

function update_value(U, I, i::Int, j::Int, k::Int, p::Int, A, env::Environment, veh::Vehicle)   
    x1 = env.x_grid[i]
    x2 = env.y_grid[j]
    x3 = env.theta_grid[k]
    x4 = env.phi_grid[p]

    x = SVector{4, Float64}(x1, x2, x3, x4)

    # compute u_ijk update
    up_min = Inf
    init_ijkp = false

    # code in this loop gets iterated the most
    for ia in 1:size(A,1)
        xdot = car_phi_EoM(x, A[ia], veh)
    
        # calculate upwind indices
        i_uw = i + Int(sign(xdot[1]))
        j_uw = j + Int(sign(xdot[2]))
        k_uw = k + Int(sign(xdot[3]))
        p_uw = p + Int(sign(xdot[4]))

        k_uw == size(env.theta_grid,1)+1 ? k_uw = 2 : k_uw = k_uw
        k_uw == 0 ? k_uw = size(env.theta_grid,1)-1 : k_uw = k_uw

        if any((I[i_uw,j,k,p], I[i,j_uw,k,p], I[i,j,k_uw,p], I[i,j,k,p_uw])) == true
            # pull value from upwind points
            u_i_uw = U[i_uw, j, k, p]
            u_j_uw = U[i, j_uw, k, p]
            u_k_uw = U[i, j, k_uw, p]
            u_p_uw = U[i, j, k, p_uw]

            # calculate value for given action
            upa_ijkp = val_eq(xdot, u_i_uw, u_j_uw, u_k_uw, u_p_uw, env)    # TO-DO: update

            if upa_ijkp < up_min
                up_min = upa_ijkp
            end

            init_ijkp = true
        end
    end

    if init_ijkp == true
        up_ijkp = up_min
    else
        up_ijkp = U[i,j,k,p]
    end

    return up_ijkp, init_ijkp
end

# TO-DO: update
function val_eq(xdot::SVector{4, Float64}, u_i_uw::Float64, u_j_uw::Float64, u_k_uw::Float64, u_p_uw::Float64, env::Environment)
    s1 = sign(xdot[1])
    s2 = sign(xdot[2])
    s3 = sign(xdot[3])
    s4 = sign(xdot[4])

    num = 1.0 + s1/env.h_xy*xdot[1]*u_i_uw + s2/env.h_xy*xdot[2]*u_j_uw + s3/env.h_theta*xdot[3]*u_k_uw + s4/env.h_phi*xdot[4]*u_p_uw
    den = s1/env.h_xy*xdot[1] + s2/env.h_xy*xdot[2] + s3/env.h_theta*xdot[3] + s4/env.h_phi*xdot[4]

    u_ijkp = num/den

    return u_ijkp
end

# main function to iteratively calculate value function using HJB
function solve_HJB_PDE(A, du_tol, max_steps, env::Environment, veh::Vehicle, anim_bool)
    # initialize U, Up, I
    U, Up, I, Ip, T, O = initialize_value_array(env, veh)

    # main function loop
    du_max = maximum(Up)
    step = 1
    gs_pass = 1
    # anim_U = @animate 
    while du_max > du_tol && step < max_steps
        sweep = gauss_seidel_sweep(gs_pass, env)
        
        for i in sweep[1]
            for j in sweep[2]
                for k in sweep[3]
                    for p in sweep[4]
                        if T[i,j,k,p] == false && O[i,j,k,p] == false
                            Up[i,j,k,p], Ip[i,j,k,p] = update_value(Up, Ip, i, j, k, p, A, env, veh)
                        end
                    end
                end
            end
        end

        Up[:,:,end,:] = deepcopy(Up[:,:,1,:])
        Ip[:,:,end,:] = deepcopy(Ip[:,:,1,:])

        # compare U and Up to check convergence
        dU = Up - U
        du_max = maximum(abs.(dU))

        U = deepcopy(Up)
        I = deepcopy(Ip)

        println("step: ", step, ", du_max = ", du_max)

        if gs_pass == 2^4
            gs_pass = 1
        else
            gs_pass += 1
        end

        step += 1

        # animation ---
        if anim_bool == true
            theta_plot = 1/2*pi
            if theta_plot in env.theta_grid
                k_plot = indexin(theta_plot, env.theta_grid)[1]
            else
                k_plot = searchsortedfirst(env.theta_grid, theta_plot) - 1
            end

            p_k = heatmap(env.x_grid, env.y_grid, transpose(U[:,:,k_plot]), clim=(0,10),
                        # xlim=(-3.5,5.5),
                        aspect_ratio=:equal, 
                        size=(750,1000),
                        # xlabel="x-axis [m]", ylabel="y-axis [m]", 
                        colorbar_title = "time-to-target [s]",
                        legend=:topright,
                        # legend=true, 
                        colorbar=false,
                        legend_font_pointsize = 11,
                        top_margin = -30*Plots.mm,
                        left_margin = 8*Plots.mm,
                        bottom_margin = 4*Plots.mm)

            plot_polygon(p_k, env.W, 3, :black, "Workspace")
            plot_polygon(p_k, env.T_xy, 3, :green, "Target Set")
            plot_polygon(p_k, env.O_vec[1], 3, :red, "Obstacle")
            for O in env.O_vec
                plot_polygon(p_k, O, 3, :red, "")
            end

            # plot vehicle figure
            x_pos = 6.75
            y_pos = 5

            x_max = x_pos + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
            y_min = y_pos - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

            x = [x_pos, y_pos, env.theta_grid[k_plot]]
            
            E_arr = pose_to_edges(x, veh)
            V = [[E_arr[1][1] E_arr[1][2]];
                [E_arr[2][1] E_arr[2][2]];
                [E_arr[3][1] E_arr[3][2]];
                [E_arr[4][1] E_arr[4][2]]]
                
            plot!(p_k, [x_max], [y_pos], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot_polygon(p_k, V, 2, :blue, "Vehicle")

            # plot step count
            annotate!(x_pos, y_pos+1.5, text("step:\n$(step-1)", 14))

            display(p_k)
        end
    end
    # gif(anim_U, algs_path*"HJB-planner/figures/hjb_growth.gif", fps=3)

    return U, T, O
end

function gauss_seidel_sweep(pass, env::Environment)
    N_x = size(env.x_grid, 1)
    N_y = size(env.y_grid, 1)
    N_theta = size(env.theta_grid, 1)
    N_phi = size(env.phi_grid, 1)

    sweep_nom = [2:1:N_x-1, 2:1:N_y-1, 1:1:N_theta-1, 2:1:N_phi-1]
    sweep_gs = deepcopy(sweep_nom)

    order = reverse(digits(pass-1, base=2, pad=4))

    for i in 1:4
        if order[i] == 1
            sweep_gs[i] = reverse(sweep_nom[i])
        end
    end

    return sweep_gs

    # # define ijk iterators for Gauss-Seidel sweeping scheme
    # gs_sweeps = [[2:N_x-1, 2:N_y-1, 1:N_theta-1],     # FFF
    #             [2:N_x-1, 2:N_y-1, reverse(1:N_theta-1)],  # FFB
    #             [2:N_x-1, reverse(2:N_y-1), 1:N_theta-1],  # FBF
    #             [reverse(2:N_x-1), 2:N_y-1, 1:N_theta-1],  # BFF
    #             [2:N_x-1, reverse(2:N_y-1), reverse(1:N_theta-1)], # FBB
    #             [reverse(2:N_x-1), reverse(2:N_y-1), 1:N_theta-1], # BBF
    #             [reverse(2:N_x-1), 2:N_y-1, reverse(1:N_theta-1)], # BFB
    #             [reverse(2:N_x-1), reverse(2:N_y-1), reverse(1:N_theta-1)]]   # BBB
end


# equations of motion for 4 DoF kinematic bicycle model
function car_phi_EoM(x::SVector{4, Float64}, u::Vector{Float64}, param::Vehicle)
    xdot1 = u[1]*cos(x[3])
    xdot2 = u[1]*sin(x[3])
    xdot3 = u[1]*(1/param.wb)*tan(x[4])
    xdot4 = u[2]

    xdot = SVector{4, Float64}(xdot1, xdot2, xdot3, xdot4)

    return xdot
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k::Vector{Float64}, u::Vector{Float64}, dt, EoM::Function, param::Vehicle)
    x_ks = SVector{4, Float64}(x_k)

    w1 = EoM(x_ks, u, param)
    w2 = EoM(x_ks + w1*dt/2, u, param)
    w3 = EoM(x_ks + w2*dt/2, u, param)
    w4 = EoM(x_ks + w3*dt, u, param)

    x_k1s = x_ks + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    # convert x_k1 back
    x_k1 = [x_k1s[1], x_k1s[2], x_k1s[3], x_k1s[4]]   # TO-DO: clean up variable types

    # adjust theta within bounds
    if x_k1[3] > pi
        x_k1[3] -= 2*pi
    elseif x_k1[3] < -pi
        x_k1[3] += 2*pi
    end

    return [x_k1[1], x_k1[2], x_k1[3], x_k1[4]]
end

function plot_polygon(my_plot, P, lw, lc, ll)
    P_x_pts = [P[:,1]; P[1,1]]
    P_y_pts = [P[:,2]; P[1,2]]

    plot!(my_plot, P_x_pts, P_y_pts, linewidth=lw, linecolor=lc, label=ll)
end



# Takei functions
# computes value update for grid point ijk using HJB finite difference scheme
function update_value_UCLA(U, i, j, k, env::Environment, veh::Vehicle)
    theta_k = env.theta_grid[k]

    xi_k = Int(sign(cos(theta_k)))
    nu_k = Int(sign(sin(theta_k)))

    ip = i + xi_k
    in = i - xi_k
    jp = j + nu_k
    jn = j - nu_k
    k == size(env.theta_grid,1)-1 ? kp = 1 : kp = k + 1
    k == 1 ? kn = size(env.theta_grid,1)-1 : kn = k - 1

    u_ip = U[ip,j,k] 
    u_in = U[in,j,k]
    u_jp = U[i,jp,k]
    u_jn = U[i,jn,k]
    u_kp = U[i,j,kp]
    u_kn = U[i,j,kn]

    # compute u_ijk update
    G_p1 = G_p1_ijk(theta_k, u_ip, u_jp, u_kp, u_kn, env, veh)
    F_p1 = F_p1_ijk(theta_k, u_ip, u_jp, env, veh)
    G_n1 = G_n1_ijk(theta_k, u_in, u_jn, u_kp, u_kn, env, veh)
    F_n1 = F_n1_ijk(theta_k, u_in, u_jn, env, veh)

    up_ijk = min(G_p1, F_p1, G_n1, F_n1, U[i,j,k])

    return up_ijk
end

# finite difference approximations
# (eq 30)
function G_p1_ijk(theta_k, u_ip, u_jp, u_kp, u_kn, env::Environment, veh::Vehicle)
    rho = veh.wb/tan(veh.c_phi)

    num = env.h_xy/veh.c_vf + abs(cos(theta_k))*u_ip + abs(sin(theta_k))*u_jp + env.h_xy/(rho*env.h_theta)*min(u_kp, u_kn)
    den = abs(cos(theta_k)) + abs(sin(theta_k)) + env.h_xy/(rho*env.h_theta)
    G_p1 = num/den

    return G_p1
end

# (eq 31)
function F_p1_ijk(theta_k, u_ip, u_jp, env::Environment, veh::Vehicle)
    num = env.h_xy/veh.c_vf + abs(cos(theta_k))*u_ip + abs(sin(theta_k))*u_jp
    den = abs(cos(theta_k)) + abs(sin(theta_k))
    F_p1 = num/den

    return F_p1
end

# (eq 32)
function G_n1_ijk(theta_k, u_in, u_jn, u_kp, u_kn, env::Environment, veh::Vehicle)
    rho = veh.wb/tan(veh.c_phi)

    num = env.h_xy/veh.c_vb + abs(cos(theta_k))*u_in + abs(sin(theta_k))*u_jn + env.h_xy/(rho*env.h_theta)*min(u_kp, u_kn)
    den = abs(cos(theta_k)) + abs(sin(theta_k)) + env.h_xy/(rho*env.h_theta)
    G_n1 = num/den

    return G_n1
end

# (eq 33)
function F_n1_ijk(theta_k, u_in, u_jn, env::Environment, veh::Vehicle)
    num = env.h_xy/veh.c_vb + abs(cos(theta_k))*u_in + abs(sin(theta_k))*u_jn
    den = abs(cos(theta_k)) + abs(sin(theta_k))
    F_n1 = num/den

    return F_n1
end