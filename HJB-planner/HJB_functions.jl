# notation:
#   x_p - horizontal position
#   y_p - vertical position
#   x - point [x_p, y_p]
#   y - pose [x, theta]
#   z - path, parameterized by time
#   sigma - speed control
#   a - steering control
#   T - target set
#   f - dynamics
#   rho - minimum turning radiuss
#   u - value function for given dynamics and target set
#   u_x, u_y, u_theta - partial derivatives of value function wrt x, y, theta

using Random
rng = MersenneTwister(1234)

# define discretized grid struct
struct Environment
    h_xy::Float64
    h_theta::Float64
    x_grid::StepRangeLen
    y_grid::StepRangeLen
    theta_grid::StepRangeLen
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

# target set checker
function in_target_set(y, env::Environment, veh::Vehicle)
    # check position of corners
    V_c = pose_to_corners(y, veh)
    rows = [collect(1:size(env.T_xy, 1)); 1]

    for c in V_c
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
    V_c = pose_to_corners(y, veh::Vehicle)

    for c in V_c 
        for Oi in env.O_vec
            rows = [collect(1:size(Oi, 1)); 1]

            ineqs = zeros(Bool, size(Oi, 1))
            for i in 1:size(Oi, 1)
                i1 = rows[i]
                i2 = rows[i+1]

                x1 = Oi[i1,1]
                y1 = Oi[i1,2]
                x2 = Oi[i2,1]
                y2 = Oi[i2,2]

                val = (y1 - y2)*c[1] + (x2 - x1)*c[2] + x1*y2 - x2*y1

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

# workspace checker
function in_workspace(y, env::Environment, veh::Vehicle)
    V_c = pose_to_corners(y, veh::Vehicle)

    for c in V_c 
        rows = [collect(1:size(env.W, 1)); 1]

        for i in 1:size(env.W, 1)
            i1 = rows[i]
            i2 = rows[i+1]

            x1 = env.W[i1,1]
            y1 = env.W[i1,2]
            x2 = env.W[i2,1]
            y2 = env.W[i2,2]

            val = (y1 - y2)*c[1] + (x2 - x1)*c[2] + x1*y2 - x2*y1

            if val < 0
                return false
            end
        end
    end

    return true
end

# calculates position of corners of rectangular vehicle
function pose_to_corners(y, veh::Vehicle)
    x_FR = y[1] + (veh.l - veh.b2a)*cos(y[3]) + (veh.w/2)*sin(y[3])
    x_FL = y[1] + (veh.l - veh.b2a)*cos(y[3]) - (veh.w/2)*sin(y[3])
    x_BR = y[1] - veh.b2a*cos(y[3]) + (veh.w/2)*sin(y[3])
    x_BL = y[1] - veh.b2a*cos(y[3]) - (veh.w/2)*sin(y[3])
    
    y_FR = y[2] + (veh.l - veh.b2a)*sin(y[3]) - (veh.w/2)*cos(y[3])
    y_FL = y[2] + (veh.l - veh.b2a)*sin(y[3]) + (veh.w/2)*cos(y[3])
    y_BR = y[2] - veh.b2a*sin(y[3]) - (veh.w/2)*cos(y[3])
    y_BL = y[2] - veh.b2a*sin(y[3]) + (veh.w/2)*cos(y[3])
    
    V_c = [[x_FR, y_FR],
            [x_FL, y_FL],
            [x_BL, y_BL],
            [x_BR, y_BR]]

    return V_c
end

# initialize value approximations
function initialize_value_array(env::Environment, veh::Vehicle)
    Up = ones(Float64, length(env.x_grid), length(env.y_grid), length(env.theta_grid))

    for i in 1:length(env.x_grid)
        for j in 1:length(env.y_grid)
            for k in 1:length(env.theta_grid)
                x_i = env.x_grid[i]
                y_j = env.y_grid[j]
                theta_k = env.theta_grid[k]

                y_ijk = [x_i, y_j, theta_k]

                if in_target_set(y_ijk, env, veh) == true
                    Up[i,j,k] = 0
                else
                    Up[i,j,k] = 50
                end
            end
        end
    end

    U = deepcopy(Up)

    return U, Up
end

# computes value update for grid point ijk using HJB finite difference scheme
function update_value(U, i, j, k, env::Environment, veh::Vehicle)
    theta_k = env.theta_grid[k]

    xi_k = Int(sign(cos(theta_k)))
    nu_k = Int(sign(sin(theta_k)))

    ip = i + xi_k
    in = i - xi_k
    jp = j + nu_k
    jn = j - nu_k
    k == length(env.theta_grid)-1 ? kp = 1 : kp = k + 1
    k == 1 ? kn = length(env.theta_grid)-1 : kn = k - 1

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

# main function to iteratively calculate value function using HJB
function solve_HJB_PDE(du_tol, max_reps, env::Environment, veh::Vehicle, animate)
    N_x = length(env.x_grid)
    N_y = length(env.y_grid)
    N_theta = length(env.theta_grid)

    # define ijk iterators for Gauss-Seidel sweeping scheme
    gs_sweeps = [[2:N_x-1, 2:N_y-1, 1:N_theta-1],     # FFF
                [2:N_x-1, 2:N_y-1, reverse(1:N_theta-1)],  # FFB
                [2:N_x-1, reverse(2:N_y-1), 1:N_theta-1],  # FBF
                [reverse(2:N_x-1), 2:N_y-1, 1:N_theta-1],  # BFF
                [2:N_x-1, reverse(2:N_y-1), reverse(1:N_theta-1)], # FBB
                [reverse(2:N_x-1), reverse(2:N_y-1), 1:N_theta-1], # BBF
                [reverse(2:N_x-1), 2:N_y-1, reverse(1:N_theta-1)], # BFB
                [reverse(2:N_x-1), reverse(2:N_y-1), reverse(1:N_theta-1)]]   # BBB

    # initialize U and Up
    (U, Up) = initialize_value_array(env, veh)

    # main function loop
    du_max = maximum(Up)
    rep = 1
    gs = 1
    # anim_U = @animate 
    while du_max > du_tol && rep < max_reps
        sweep = gs_sweeps[gs]
        
        for i in sweep[1]
            for j in sweep[2]
                for k in sweep[3]
                    x_i = env.x_grid[i]
                    y_j = env.y_grid[j]
                    theta_k = env.theta_grid[k]

                    y_ijk = [x_i, y_j, theta_k]
                    
                    if in_obstacle_set(y_ijk, env, veh) == false
                        Up[i,j,k] = update_value(U, i, j, k, env, veh)
                    end
                end
            end
        end

        Up[:,:,end] = deepcopy(Up[:,:,1])

        # compare U and Up to check convergence
        dU = Up - U
        du_max = maximum(abs.(dU))
        println(du_max)

        U = deepcopy(Up)

        if gs == 8
            gs = 1
        else
            gs += 1
        end

        rep += 1

        # animation ---
        if animate == true
            theta_plot = pi/2
            if theta_plot in env.theta_grid
                k_plot = indexin(theta_plot, env.theta_grid)[1]
            else
                k_plot = searchsortedfirst(env.theta_grid, theta_plot) - 1
            end

            p_k = heatmap(env.x_grid, env.y_grid, transpose(U[:,:,k_plot]), clim=(0,15),
                        aspect_ratio=:equal, size=(1000,850),
                        xlabel="x-axis [m]", ylabel="y-axis [m]", 
                        title="HJB Value Function",
                        titlefontsize = 20,
                        colorbar_title = "time-to-target [s]",
                        legend=:topright,
                        legend_font_pointsize = 11,
                        top_margin = -30*Plots.mm,
                        left_margin = 8*Plots.mm,
                        bottom_margin = -8*Plots.mm)

            plot_polygon(p_k, env.W, 2, :black, "Workspace")
            plot_polygon(p_k, env.T_xy, 2, :green, "Target Set")
            plot_polygon(p_k, env.O_vec[1], 2, :red, "Obstacle")    # TO-DO: make obstacle plotting cleaner
            plot_polygon(p_k, env.O_vec[2], 2, :red, "")
            plot_polygon(p_k, env.O_vec[3], 2, :red, "")

            # plot vehicle figure
            y = [6, -4, env.theta_grid[k_plot]]
            
            V_c = pose_to_corners(y, unit_car)
            V = [[V_c[1][1] V_c[1][2]];
                [V_c[2][1] V_c[2][2]];
                [V_c[3][1] V_c[3][2]];
                [V_c[4][1] V_c[4][2]]]
                
            plot!(p_k, [x_max], [-4], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot!(p_k, [6], [-4], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot_polygon(p_k, V, 2, :blue, "Vehicle Orientation")

            # plot step count
            annotate!(6, 0, text("step:\n$(rep-1)", 14))
        end
    end
    # gif(anim_U, "hjb_90_growth.gif", fps=6)

    return U
end

# interpolate U to return a value for any point (x,y,theta)
function interp_value(y, U, env::Environment)
    # adjust theta within bounds
    if y[3] > pi
        y[3] -= 2*pi
    elseif y[3] < -pi
        y[3] += 2*pi
    end

    # implements trilinear interpolation from Wikipedia
    println(y[1])
    if y[1] in env.x_grid
        i_0 = indexin(y[1], env.x_grid)[1]
    else
        i_0 = searchsortedfirst(env.x_grid, y[1]) - 1
    end

    if y[2] in env.y_grid
        j_0 = indexin(y[2], env.y_grid)[1]
    else
        j_0 = searchsortedfirst(env.y_grid, y[2]) - 1
    end

    if y[3] in env.theta_grid
        k_0 = indexin(y[3], env.theta_grid)[1]
    else
        k_0 = searchsortedfirst(env.theta_grid, y[3]) - 1
    end

    i_1 = i_0 + 1
    j_1 = j_0 + 1
    k_0 == length(env.theta_grid) ? k_1 = 1 : k_1 = k_0 + 1

    x_0 = env.x_grid[i_0]           # ISSUE: indexing to 0 (i_0 assigned in line 379)
    y_0 = env.y_grid[j_0]           # ISSUE: indexing to 0 (j_0 assigned in line 385)
    theta_0 = env.theta_grid[k_0]   # ISSUE: indexing to 0

    x_1 = env.x_grid[i_1]
    y_1 = env.y_grid[j_1]
    theta_1 = env.theta_grid[k_1]

    x_d = (y[1] - x_0)/(x_1 - x_0)
    y_d = (y[2] - y_0)/(y_1 - y_0)
    theta_d = (y[3] - theta_0)/(theta_1 - theta_0)

    u_000 = U[i_0, j_0, k_0]
    u_100 = U[i_1, j_0, k_0]
    u_010 = U[i_0, j_1, k_0]
    u_110 = U[i_1, j_1, k_0]
    u_001 = U[i_0, j_0, k_1]
    u_101 = U[i_1, j_0, k_1]
    u_011 = U[i_0, j_1, k_1]
    u_111 = U[i_1, j_1, k_1]

    u_00 = u_000*(1 - x_d) + u_100*x_d
    u_01 = u_001*(1 - x_d) + u_101*x_d
    u_10 = u_010*(1 - x_d) + u_110*x_d
    u_11 = u_011*(1 - x_d) + u_111*x_d

    u_0 = u_00*(1 - y_d) + u_10*y_d
    u_1 = u_01*(1 - y_d) + u_11*y_d

    u_itp = u_0*(1 - theta_d) + u_1*theta_d

    return u_itp
end

function HJB_planner(y_0, U, dt, env::Environment, veh::Vehicle)
    max_steps = 5000

    if typeof(y_0) != Vector{Float64}
        y_0 = convert(Vector{Float64}, y_0)
    end

    y_k = y_0
    y_path = [y_k]  

    u_path = []

    step = 0
    while in_target_set(y_k, env, veh) == false && step < max_steps
        step += 1

        # calculate optimal action
        u_k = optimal_action_HJB(y_k, U, env, veh)
        push!(u_path, u_k)

        # set noise parameters
        std_v = 0
        std_phi = 0

        # simulate forward one time step
        y_k1 = runge_kutta_4(car_EoM, y_k, u_k, std_v, std_phi, dt, veh)
        push!(y_path, y_k1)

        y_k = deepcopy(y_k1)
    end

    println("steps in HJB path: ", step)

    return y_path, u_path, step
end



# NOTE: want to understand why value iteration gets hung up on certain values

# For Himanshu:
    # TO-DO: remove backwards motion
    # TO-DO: time the optimal_action
    # TO-DO:

    # return to POMDP: optimal action at given state

    # dt = 0.1 vs dt = 1.0, what issues will appear? hitting obstacles? 
    # reactive controller will modify velocity on path, how will this affect path execution and safety?

    # delta rollout action

    # calculate limited horizon path ahead of car, based on current velocity
    #   - generate next 30 points and actions

# calculates optimal action as a function of state
function optimal_action_HJB(y, U, env::Environment, veh::Vehicle)
    # adjust position when near edge for boundary issues
    if abs(y[1] - env.x_grid[1]) < 2*env.h_xy
        y[1] = env.x_grid[3]
    elseif abs(y[1] - env.x_grid[end]) < 2*env.h_xy
        y[1] = env.x_grid[end-2]
    end

    if abs(y[2] - env.y_grid[1]) < 2*env.h_xy
        y[2] = env.y_grid[3]
    elseif abs(y[2] - env.y_grid[end]) < 2*env.h_xy
        y[2] = env.y_grid[end-2]
    end

    # approximate partial derivatives using centered difference scheme
    du_dx = (interp_value(y+[env.h_xy,0,0], U, env) - interp_value(y-[env.h_xy,0,0], U, env))/(2*env.h_xy)
    du_dy = (interp_value(y+[0,env.h_xy,0], U, env) - interp_value(y-[0,env.h_xy,0], U, env))/(2*env.h_xy)
    du_dtheta = (interp_value(y+[0,0,env.h_theta], U, env) - interp_value(y-[0,0,env.h_theta], U, env))/(2*env.h_theta)

    # calculate optimal action
    A = [[veh.c_vf, 0],
        [veh.c_vf, veh.c_phi],
        [veh.c_vf, -veh.c_phi],
        [-veh.c_vb, 0],
        [-veh.c_vb, veh.c_phi],
        [-veh.c_vb, -veh.c_phi]]

    # TO-DO: need to bias straight actions (set some margin)
    HJB_min(a) = a[1]*(du_dx*cos(y[3]) + du_dy*sin(y[3]) + du_dtheta*(1/veh.wb)*tan(a[2]))
    a_opt = argmin(HJB_min, A)

    return a_opt
end

# equations of motion for 3 DoF kinematic bicycle model
function car_EoM(x, u, std_v, std_phi, param)
    w_v = std_v*randn(rng, Float64) .+ 0
    w_phi = std_phi*randn(rng, Float64) .+ 0

    x_dot = [(u[1]+w_v)*cos(x[3]),
            (u[1]+w_v)*sin(x[3]),
            (u[1]+w_v)*(1/param.wb)*tan(u[2]+w_phi)]

    return x_dot
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(EoM::Function, x_k, u, std_v, std_phi, dt, param)
    w1 = EoM(x_k, u, std_v, std_phi, param)
    w2 = EoM(x_k + w1*dt/2, u, std_v, std_phi, param)
    w3 = EoM(x_k + w2*dt/2, u, std_v, std_phi, param)
    w4 = EoM(x_k + w3*dt, u, std_v, std_phi, param)

    x_k1 = x_k + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    # adjust theta within bounds
    if x_k1[3] > pi
        x_k1[3] -= 2*pi
    elseif x_k1[3] < -pi
        x_k1[3] += 2*pi
    end

    return x_k1
end

function plot_polygon(my_plot, P, lw, lc, ll)
    P_x_pts = [P[:,1]; P[1,1]]
    P_y_pts = [P[:,2]; P[1,2]]

    plot!(my_plot, P_x_pts, P_y_pts, linewidth=lw, linecolor=lc, label=ll)
end