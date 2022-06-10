using LinearAlgebra
using Random
rng = MersenneTwister(1234)

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"

# define discretized grid struct
struct Environment
    h_xy::Float64
    h_theta::Float64
    h_v::Float64
    x_grid::StepRangeLen
    y_grid::StepRangeLen
    theta_grid::StepRangeLen
    v_grid::StepRangeLen
    W::Matrix{Float64}
    T_xy::Matrix{Float64}
    T_theta::Vector{Vector{Float64}}
    T_v::Vector{Float64}
    O_vec::Vector{Matrix{Float64}}
end

struct Vehicle
    c_vf::Float64   # max forward speed [m/s]
    c_vb::Float64   # max backward speed [m/s]
    c_af::Float64   # max forward accel [m/s^2]
    c_ab::Float64   # max backward accel [m/s^2]
    c_phi::Float64  # max steering angle [rad]
    wb::Float64     # wheelbase [m]
    l::Float64      # length [m]
    w::Float64      # width [m]
    b2a::Float64    # rear bumber to rear axle [m]
end

# ?: should velocity limits be used to define state space bounds?
#   - seems like it could be more of an "obstacle" in the state space, rather than the defined edge
#   - consider that x,y(,theta?) are defined by environment, v is defined by vehicle
#   - BUT: axes of state space are not physical representation of environment, they represent the state of the vehicle
#       - so think makes sense as defined boundary

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

    # check velocity
    if y[4] < minimum(env.T_v) || y[4] > maximum(env.T_v)
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

# state space checker
# ?: is this actually used somewhere?
function in_state_space(y, env::Environment, veh::Vehicle)
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

            if val <= 0
                return false
            end
        end
    end

    # check velocity
    if y[4] < minimum(env.T_v) || y[4] > maximum(env.T_v)
        return false
    end

    return true
end

# ISSUE: want boundary velocities to still be valid
#   - want to make boundaries inclusive (including for xy points)
function on_boundary(i, j, k, p, env::Environment)
    if (i == 1 || i == size(env.x_grid,1)) || (j == 1 || j == size(env.y_grid,1)) || (p == 1 || p == size(env.v_grid,1))
        return true
    else
        return false
    end
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
    Up = zeros(Float64, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.v_grid,1)))
    Ip = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.v_grid,1)))
    T = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.v_grid,1)))
    O = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1), size(env.v_grid,1)))

    for i in 1:size(env.x_grid,1)
        for j in 1:size(env.y_grid,1)
            for k in 1:size(env.theta_grid,1)
                for p in 1:size(env.v_grid,1)
                    x_i = env.x_grid[i]
                    y_j = env.y_grid[j]
                    theta_k = env.theta_grid[k]
                    v_p = env.v_grid[p]

                    y_ijk = [x_i, y_j, theta_k, v_p]

                    if in_target_set(y_ijk, env, veh) == true
                        Up[i,j,k,p] = 0.0
                        Ip[i,j,k,p] = true
                        T[i,j,k,p] = true

                    # ISSUE: initializes everything at 30.0 -> all on boundary?
                    elseif on_boundary(i, j, k, p, env) == true #|| in_obstacle_set(y_ijk, env, veh) == true
                        Up[i,j,k,p] = 30.0
                        Ip[i,j,k,p] = false
                        O[i,j,k,p] = true
                    else
                        # Up[i,j,k] = norm([env.x_grid[i] - 0, env.y_grid[j] - 0])
                        Up[i,j,k,p] = 20.0
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

# ISSUE: looks like memory allocation (285 occurences) is biggest difference from Takei
#   - most were 285 in first step, believe these were non-initialized states
#   - a few were 315, possibly states with new info (entered if-true statement)
#   - actually more allocations without x= (285 -> 300) ???
#   - issue to recreate up_arr, etc every time function runs? (new allocation?)

#   - down to 283,297 allocations
#   - down to 115,127 allocations

#   - spikes to  290 allocations when node is initialized (so, most of the time) (or only 139?)
#   - runtimes get longer over time, even with same number of allocations

#   - Static Arrays?

# a = @allocated begin
#     Block to test
# end; if a > 0 println(a) end

#   - look at type stability
function update_value(U, I, i::Int, j::Int, k::Int, p::Int, A, EoM::Function, env::Environment, veh::Vehicle)
    theta_k = env.theta_grid[k]

    x = [env.x_grid[i],
        env.y_grid[j],
        theta_k,
        env.v_grid[p]]

    # compute u_ijk update
    up_arr = zeros(Float64, size(A,1))
    init_ijkp = false

    # ISSUE: 60 allocations in for loop
    #   - 22 allocations for each iteration
    #   - only needs to keep track of up_arr, init_ijk
    # @time begin
    for ia in 1:size(A,1)
        a_a = A[ia][1]
        a_phi = A[ia][2]

        xdot = EoM(x, [a_a, a_phi], veh)
    
        # calculate upwind indices
        i_uw = i + Int(sign(xdot[1]))   # 2 allocations
        j_uw = j + Int(sign(xdot[2]))
        k_uw = k + Int(sign(xdot[3]))
        p_uw = p + Int(sign(xdot[4]))

        k_uw == size(env.theta_grid,1)+1 ? k_uw = 2 : k_uw = k_uw
        k_uw == 0 ? k_uw = size(env.theta_grid,1)-1 : k_uw = k_uw

        if any([I[i_uw,j,k,p], I[i,j_uw,k,p], I[i,j,k_uw,p], I[i,j,k,p_uw]]) == true
            # @time begin # 4 allocations
            # pull value from upwind points
            u_i_uw = U[i_uw, j, k, p]
            u_j_uw = U[i, j_uw, k, p]
            u_k_uw = U[i, j, k_uw, p]
            u_p_uw = U[i, j, k, p_uw]

            # calculate value for given action
            up_ijkp = val_eq(xdot, u_i_uw, u_j_uw, u_k_uw, u_p_uw, env)  # 1 allocation

            up_arr[ia] = up_ijkp
            init_ijkp = true
            # end
        else
            up_arr[ia] = Inf
            continue
        end
    # end
    end

    if init_ijkp == true
        up_ijkp = minimum(up_arr)
    else
        up_ijkp = U[i,j,k,p]
    end

    return up_ijkp, init_ijkp
end

function val_eq(xdot::Vector{Float64}, u_i_uw::Float64, u_j_uw::Float64, u_k_uw::Float64, u_p_uw::Float64, env::Environment)
    s1 = sign(xdot[1])
    s2 = sign(xdot[2])
    s3 = sign(xdot[3])
    s4 = sign(xdot[4])

    num = 1.0 + s1/env.h_xy*xdot[1]*u_i_uw + s2/env.h_xy*xdot[2]*u_j_uw + s3/env.h_theta*xdot[3]*u_k_uw + s4/env.h_v*xdot[4]*u_p_uw
    den = s1/env.h_xy*xdot[1] + s2/env.h_xy*xdot[2] + s3/env.h_theta*xdot[3] + s4/env.h_v*xdot[4]

    u_ijkp = num/den

    return u_ijkp
end

# main function to iteratively calculate value function using HJB
function solve_HJB_PDE(du_tol, max_steps, env::Environment, veh::Vehicle, anim_bool)
    N_x = size(env.x_grid,1)
    N_y = size(env.y_grid,1)
    N_theta = size(env.theta_grid,1)
    N_v = size(env.v_grid,1)

    # define ijk iterators for Gauss-Seidel sweeping scheme
    gs_sweeps = [[2:N_x-1, 2:N_y-1, 1:N_theta-1, 2:N_v-1],     # FFFF
                [2:N_x-1, 2:N_y-1, reverse(1:N_theta-1), 2:N_v-1],  # FFBF
                [2:N_x-1, reverse(2:N_y-1), 1:N_theta-1, 2:N_v-1],  # FBFF
                [reverse(2:N_x-1), 2:N_y-1, 1:N_theta-1, 2:N_v-1],  # BFFF
                [2:N_x-1, reverse(2:N_y-1), reverse(1:N_theta-1), 2:N_v-1], # FBBF
                [reverse(2:N_x-1), reverse(2:N_y-1), 1:N_theta-1, 2:N_v-1], # BBFF
                [reverse(2:N_x-1), 2:N_y-1, reverse(1:N_theta-1), 2:N_v-1], # BFBF
                [reverse(2:N_x-1), reverse(2:N_y-1), reverse(1:N_theta-1), 2:N_v-1],   # BBBF
                [2:N_x-1, 2:N_y-1, 1:N_theta-1, reverse(2:N_v-1)],     # FFFB
                [2:N_x-1, 2:N_y-1, reverse(1:N_theta-1), reverse(2:N_v-1)],  # FFBB
                [2:N_x-1, reverse(2:N_y-1), 1:N_theta-1, reverse(2:N_v-1)],  # FBFB
                [reverse(2:N_x-1), 2:N_y-1, 1:N_theta-1, reverse(2:N_v-1)],  # BFFB
                [2:N_x-1, reverse(2:N_y-1), reverse(1:N_theta-1), reverse(2:N_v-1)], # FBBB
                [reverse(2:N_x-1), reverse(2:N_y-1), 1:N_theta-1, reverse(2:N_v-1)], # BBFB
                [reverse(2:N_x-1), 2:N_y-1, reverse(1:N_theta-1), reverse(2:N_v-1)], # BFBB
                [reverse(2:N_x-1), reverse(2:N_y-1), reverse(1:N_theta-1), reverse(2:N_v-1)]] # BBBB

    A = [[a_a,a_phi] for a_a in [-veh.c_ab, 0.0, veh.c_af], a_phi in [-veh.c_phi, 0.0, veh.c_phi]]
    A = reshape(A, (9,1))

    # initialize U, Up, I
    U, Up, I, Ip, T, O = initialize_value_array(env, veh)

    # for (p,v) in enumerate(env.v_grid)
    #     println("v = ", v)
    #     display(U[:,:,1,p])
    # end

    # test_ids = [[2, 2, 16],
    #             [3, 3, 16],
    #             [4, 4, 16],
    #             [5, 5, 16],
    #             [6, 6, 16],
    #             [7, 7, 16],
    #             [8, 8, 16]]

    test_ids = [[9, 2, 28],
                [9, 3, 28],
                [9, 4, 28],
                [9, 5, 28],
                [9, 6, 28],
                [9, 7, 28],
                [9, 8, 28]]

    # test_ids = [[2, 2, 16],
    #             [6, 3, 21],
    #             [12, 8, 6],
    #             [5, 5, 16]]

    # U_hist = []
    # row = []
    # for id in test_ids
    #     U_id = U[id[1], id[2], id[3]]
    #     push!(row, U_id)
    # end
    # push!(U_hist, row)
    # row = []

    # main function loop
    du_max = maximum(Up)
    step = 1
    gs = 1
    # anim_U = @animate 
    while du_max > du_tol && step < max_steps
        sweep = gs_sweeps[gs]
        
        for i in sweep[1]
            for j in sweep[2]
                for k in sweep[3]
                    for p in sweep[4]
                        if T[i,j,k,p] == false && O[i,j,k,p] == false    
                            Up[i,j,k,p], Ip[i,j,k,p] = update_value(U, I, i, j, k, p, A, car4_EoM, env, veh)
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
        
        # for id in test_ids
        #     U_id = U[id[1], id[2], id[3]]
        #     push!(row, U_id)
        # end
        # push!(U_hist, row)
        # row = []

        println("step: ", step, ", du_max = ", du_max)

        if gs == 2^4
            gs = 1
        else
            gs += 1
        end

        step += 1

        # animation ---
        if anim_bool == true
            theta_plot = pi/2
            if theta_plot in env.theta_grid
                k_plot = indexin(theta_plot, env.theta_grid)[1]
            else
                k_plot = searchsortedfirst(env.theta_grid, theta_plot) - 1
            end

            v_plot = 0.75
            if v_plot in env.v_grid
                p_plot = indexin(v_plot, env.v_grid)[1]
            else
                p_plot = searchsortedfirst(env.v_grid, v_plot) - 1
            end

            p_k = heatmap(env.x_grid, env.y_grid, transpose(U[:,:,k_plot,p_plot]), clim=(0,30),
                        # xlim=(-3.5,5.5),
                        aspect_ratio=:equal, 
                        # size=(775,1050),
                        size=(800,1050),
                        xlabel="x-axis [m]", ylabel="y-axis [m]", 
                        title="HJB Value Function",
                        titlefontsize = 20,
                        colorbar_title = "time-to-target [s]",
                        # legend=:topright,
                        legend=false, colorbar=false,
                        legend_font_pointsize = 11,
                        top_margin = -30*Plots.mm,
                        left_margin = 8*Plots.mm,
                        bottom_margin = -8*Plots.mm)

            # for i in 1:size(env.x_grid,1)
            #     for j in 1:size(env.y_grid,1)
            #         if I[i,j,k_plot] == true
            #             plot!(p_k, [env.x_grid[i]], [env.y_grid[j]],
            #                 markershape=:circle, markercolor=:white, markersize=3, 
            #                 label = "")
            #         end
            #     end
            # end
            plot!(p_k)

            plot_polygon(p_k, env.W, 2, :black, "Workspace")
            plot_polygon(p_k, env.T_xy, 2, :green, "Target Set")
            # plot_polygon(p_k, env.O_vec[1], 2, :red, "Obstacle")    # TO-DO: make obstacle plotting cleaner
            # plot_polygon(p_k, env.O_vec[2], 2, :red, "")
            # plot_polygon(p_k, env.O_vec[3], 2, :red, "")

            # plot vehicle figure
            x_pos = 4.5
            y_pos = -1

            x_max = x_pos + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
            y_min = y_pos - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

            y = [x_pos, y_pos, env.theta_grid[k_plot], env.v_grid[p_plot]]
            
            V_c = pose_to_corners(y, unit_car)
            V = [[V_c[1][1] V_c[1][2]];
                [V_c[2][1] V_c[2][2]];
                [V_c[3][1] V_c[3][2]];
                [V_c[4][1] V_c[4][2]]]
                
            plot!(p_k, [x_max], [y_pos], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot_polygon(p_k, V, 2, :blue, "Vehicle Orientation")

            # plot step count
            annotate!(x_pos, y_pos+2, text("step:\n$(step-1)", 14))

            display(p_k)
        end
    end
    # gif(anim_U, algs_path*"HJB-planner/figures/hjb_growth.gif", fps=6)

    # p_conv = plot(xlabel="step", ylabel="value", title="HJB Value Convergence vs Step")
    # for i in 1:size(test_ids,1)
    #     id = test_ids[i]
    #     println("id: ", id)
    #     plot!(p_conv, 0:step-1, getindex.(U_hist,i), label="id = $id")
    # end
    # display(p_conv)

    display(U)

    return U
end

# may produce asymmetry due to interpolation index selection
function interp_value(y, U, env::Environment)
    y_itp = deepcopy(y)

    # map out-of-bounds xy states to boundary
    if y_itp[1] > env.x_grid[end-1]
        y_itp[1] = env.x_grid[end-1]
    elseif y_itp[1] < env.x_grid[2]
        y_itp[1] = env.x_grid[2]
    end

    if y_itp[2] > env.y_grid[end-1]
        y_itp[2] = env.y_grid[end-1]
    elseif y_itp[2] < env.y_grid[2]
        y_itp[2] = env.y_grid[2]
    end
    
    # adjust theta within bounds
    if y_itp[3] > pi
        y_itp[3] -= 2*pi
    elseif y_itp[3] < -pi
        y_itp[3] += 2*pi
    end

    # implements trilinear interpolation from Wikipedia
    # println(y[1])
    if y_itp[1] in env.x_grid
        i_0 = indexin(y_itp[1], env.x_grid)[1]
    else
        i_0 = searchsortedfirst(env.x_grid, y_itp[1]) - 1
    end

    if y_itp[2] in env.y_grid
        j_0 = indexin(y_itp[2], env.y_grid)[1]
    else
        j_0 = searchsortedfirst(env.y_grid, y_itp[2]) - 1
    end

    if y_itp[3] in env.theta_grid
        k_0 = indexin(y_itp[3], env.theta_grid)[1]
    else
        k_0 = searchsortedfirst(env.theta_grid, y_itp[3]) - 1
    end

    i_1 = i_0 + 1
    j_1 = j_0 + 1
    k_0 == size(env.theta_grid,1) ? k_1 = 1 : k_1 = k_0 + 1

    x_0 = env.x_grid[i_0]
    y_0 = env.y_grid[j_0]
    theta_0 = env.theta_grid[k_0]

    x_1 = env.x_grid[i_1]   
    y_1 = env.y_grid[j_1]
    theta_1 = env.theta_grid[k_1]

    x_d = (y_itp[1] - x_0)/(x_1 - x_0)
    y_d = (y_itp[2] - y_0)/(y_1 - y_0)
    theta_d = (y_itp[3] - theta_0)/(theta_1 - theta_0)

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

        # simulate forward one time step
        y_k1 = runge_kutta_4(y_k, u_k, dt, car4_EoM, veh)    # needs K_sub
        push!(y_path, y_k1)

        y_k = deepcopy(y_k1)
    end

    # println("steps in HJB path: ", step)

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
function car4_EoM(x, u, param::Vehicle)
    xdot = [x[4]*cos(x[3]),
            x[4]*sin(x[3]),
            x[4]*(1/param.wb)*tan(u[2]),
            u[1]]

    return xdot
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k, u, dt, EoM::Function, param::Vehicle)
    w1 = EoM(x_k, u, param)
    w2 = EoM(x_k + w1*dt/2, u, param)
    w3 = EoM(x_k + w2*dt/2, u, param)
    w4 = EoM(x_k + w3*dt, u, param)

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