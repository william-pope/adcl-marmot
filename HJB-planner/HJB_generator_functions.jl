# HJB_generator_functions.jl

using StaticArrays
using LinearAlgebra
using Rotations
using BenchmarkTools
using ProfileView

include("HJB_utils.jl")
include("dynamics_models.jl")

algs_path_mac = "/value_arraysers/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"

# main function to iteratively calculate conotinuous value function for HJB using finite difference method (STLC models only)
function solve_HJB_PDE(actions, dval_tol, max_steps, env::Environment, veh::Vehicle, EoM::Function, anim_bool)
    N_x = size(env.x_grid,1)
    N_y = size(env.y_grid,1)
    N_theta = size(env.theta_grid,1)

    # define ijk iterators for Gauss-Seidel sweeping scheme
    gs_sweeps = [[2:N_x-1, 2:N_y-1, 1:N_theta-1],     # FFF
                [2:N_x-1, 2:N_y-1, reverse(1:N_theta-1)],  # FFB
                [2:N_x-1, reverse(2:N_y-1), 1:N_theta-1],  # FBF
                [reverse(2:N_x-1), 2:N_y-1, 1:N_theta-1],  # BFF
                [2:N_x-1, reverse(2:N_y-1), reverse(1:N_theta-1)], # FBB
                [reverse(2:N_x-1), reverse(2:N_y-1), 1:N_theta-1], # BBF
                [reverse(2:N_x-1), 2:N_y-1, reverse(1:N_theta-1)], # BFB
                [reverse(2:N_x-1), reverse(2:N_y-1), reverse(1:N_theta-1)]]   # BBB

    # initialize value_array, value_array, init_array
    value_array, init_array, target_array, obstacle_array = initialize_value_array(env, veh)
    value_array_old = deepcopy(value_array)

    # main function loop
    dval_max = maximum(value_array)
    step = 1
    gs = 1
    # anim_value_array = @animate 
    while dval_max > dval_tol && step < max_steps
        sweep = gs_sweeps[gs]
        
        for i in sweep[1]
            for j in sweep[2]
                for k in sweep[3]
                    if target_array[i,j,k] == false && obstacle_array[i,j,k] == false
                        # # finite difference method
                        # value_array[i,j,k], init_array[i,j,k] = update_node_value_FDM(value_array, init_array, i, j, k, actions, env, veh, EoM)
                        
                        # semi-Lagrangian method
                        value_array[i,j,k] = update_node_value_SL(value_array, i, j, k, actions, env, veh, EoM)
                    end
                end
            end
        end

        value_array[:,:,end] = deepcopy(value_array[:,:,1])
        init_array[:,:,end] = deepcopy(init_array[:,:,1])

        # compare value_array and value_array to check convergence
        dval = value_array - value_array_old
        dval_max = maximum(abs.(dval))

        value_array_old = deepcopy(value_array)

        println("step: ", step, ", dval_max = ", dval_max)

        if gs == 8
            gs = 1
        else
            gs += 1
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

            p_k = heatmap(env.x_grid, env.y_grid, transpose(value_array[:,:,k_plot]), clim=(0,15),
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
            for O_vec in env.O_vec
                plot_polygon(p_k, O_vec, 3, :red, "")
            end

            # plot vehicle figure
            x_pos = 9.0
            y_pos = 4.5

            x_max = x_pos + sqrt((veh.axle_l-veh.ext2axle)^2 + (veh.ext_w/2)^2)
            y_min = y_pos - sqrt((veh.axle_l-veh.ext2axle)^2 + (veh.ext_w/2)^2)

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
    # gif(anim_value_array, algs_path*"HJB-planner/figures/hjb_growth.gif", fps=3)

    return value_array, target_array, obstacle_array
end

# initialize value approximations
function initialize_value_array(env::Environment, veh::Vehicle)
    value_array = zeros(Float64, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1)))
    init_array = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1)))
    target_array = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1)))
    obstacle_array = zeros(Bool, (size(env.x_grid,1), size(env.y_grid,1), size(env.theta_grid,1)))

    for i in 1:size(env.x_grid,1)
        for j in 1:size(env.y_grid,1)
            for k in 1:size(env.theta_grid,1)
                x_i = env.x_grid[i]
                y_j = env.y_grid[j]
                theta_k = env.theta_grid[k]

                y_ijk = [x_i, y_j, theta_k]

                if in_target_set(y_ijk, env, veh) == true
                    value_array[i,j,k] = 0.0
                    init_array[i,j,k] = true
                    target_array[i,j,k] = true
                elseif on_boundary(i, j, k, env) == true || in_obstacle_set(y_ijk, env, veh) == true
                    value_array[i,j,k] = 1000.0
                    init_array[i,j,k] = false
                    obstacle_array[i,j,k] = true
                else
                    value_array[i,j,k] = 100.0
                    init_array[i,j,k] = false
                    target_array[i,j,k] = false
                    obstacle_array[i,j,k] = false
                end
            end
        end
    end

    value_array = deepcopy(value_array)
    init_array = deepcopy(init_array)

    return value_array, init_array, target_array, obstacle_array
end

# use modules
# Julia workflows

# @code_warntype

# TO-DO: may need to do collision checking along propagation path to make sure branch doesn't jump over obstacle

# NOTE: may need to make changes to init_array process
function update_node_value_SL(value_array, i::Int, j::Int, k::Int, actions, env::Environment, veh::Vehicle, EoM::Function)
    rho = veh.axle_l/tan(veh.u_phi_max)
    
    # ISSUE: value array behaves strangely when dt is small (~<= 0.1) (seems unexpected)

    dt = 0.4    # needs to be function of grid size and vehicle speed

    x_ijk = [env.x_grid[i], env.y_grid[j], env.theta_grid[k]]

    val_p_min = Inf
    for u in actions
        x_p = runge_kutta_4(x_ijk, u, dt, EoM, veh)
        val_p = interp_value_SL(x_p, value_array, env)

        if val_p < val_p_min
            val_p_min = val_p
        end
    end

    val_ijk = dt + val_p_min

    return val_ijk
end

# TO-DO: ignoring obstacle issues for now, will need to address
#   - may need to keep track of RIC (kinda like H-J reachability stuff...)
function interp_value_SL(x::Vector{Float64}, value_array::Array{Float64, 3}, env::Environment)
    x_itp = deepcopy(x)

    # ISSUE: interpolating near boundary has same problem as near obstacles
    #   - obstacle value should not be treated as valid (but it is...)

    # adjust x,y within bounds
    x_itp[1] <= env.x_grid[1] ? x_itp[1] = env.x_grid[1]+1/16*env.h_xy : x_itp[1] = x_itp[1]
    x_itp[1] >= env.x_grid[end] ? x_itp[1] = env.x_grid[end]-1/16*env.h_xy : x_itp[1] = x_itp[1]

    x_itp[2] <= env.y_grid[1] ? x_itp[2] = env.y_grid[1]+1/16*env.h_xy : x_itp[2] = x_itp[2]
    x_itp[2] >= env.y_grid[end] ? x_itp[2] = env.y_grid[end]-1/16*env.h_xy : x_itp[2] = x_itp[2]

    # if x_itp[1] <= env.x_grid[1] || x_itp[1] >= env.x_grid[end] || x_itp[2] <= env.y_grid[1] || x_itp[2] >= env.y_grid[end]
    #     val_itp = value_array[1,1,1]
    #     return 
    # end

    # get box node indices
    i_0 = find_idx(x_itp[1], env.x_grid)
    j_0 = find_idx(x_itp[2], env.y_grid)
    k_0 = find_idx(x_itp[3], env.theta_grid)

    # NOTE: might never interp value from exactly on line, since all interp states come from propagation
    #   - would be able to simplify find_idx()

    i_1 = i_0 + 1
    j_1 = j_0 + 1
    k_0 == size(env.theta_grid,1) ? k_1 = 1 : k_1 = k_0 + 1

    # println("")
    # println("x = ", x)
    # println("x_itp = ", x_itp)
    # println("[i_0, j_0, k_0] = ", [i_0, j_0, k_0])

    # get box node states
    x_0 = env.x_grid[i_0]
    y_0 = env.y_grid[j_0]
    theta_0 = env.theta_grid[k_0]

    x_1 = env.x_grid[i_1]   
    y_1 = env.y_grid[j_1]
    theta_1 = env.theta_grid[k_1]

    x_d = (x_itp[1] - x_0)/(x_1 - x_0)
    y_d = (x_itp[2] - y_0)/(y_1 - y_0)
    theta_d = (x_itp[3] - theta_0)/(theta_1 - theta_0)

    val_000 = value_array[i_0, j_0, k_0]
    val_100 = value_array[i_1, j_0, k_0]
    val_010 = value_array[i_0, j_1, k_0]
    val_110 = value_array[i_1, j_1, k_0]
    val_001 = value_array[i_0, j_0, k_1]
    val_101 = value_array[i_1, j_0, k_1]
    val_011 = value_array[i_0, j_1, k_1]
    val_111 = value_array[i_1, j_1, k_1]

    val_00 = val_000*(1 - x_d) + val_100*x_d
    val_01 = val_001*(1 - x_d) + val_101*x_d
    val_10 = val_010*(1 - x_d) + val_110*x_d
    val_11 = val_011*(1 - x_d) + val_111*x_d

    val_0 = val_00*(1 - y_d) + val_10*y_d
    val_1 = val_01*(1 - y_d) + val_11*y_d

    val_itp = val_0*(1 - theta_d) + val_1*theta_d

    return val_itp
end

function update_node_value_FDM(value_array, init_array, i::Int, j::Int, k::Int, actions, env::Environment, veh::Vehicle, EoM::Function)   
    x1 = env.x_grid[i]
    x2 = env.y_grid[j]
    x3 = env.theta_grid[k]

    x = SVector{3, Float64}(x1, x2, x3)

    # compute val_ijk update
    val_ijk_min = Inf
    init_ijk = false

    # code in this loop gets iterated the most
    for u in actions
        xdot = EoM(x, u, veh)
    
        # calculate upwind indices
        i_uw = i + Int(sign(xdot[1]))
        j_uw = j + Int(sign(xdot[2]))
        k_uw = k + Int(sign(xdot[3]))

        k_uw == size(env.theta_grid,1)+1 ? k_uw = 2 : k_uw = k_uw
        k_uw == 0 ? k_uw = size(env.theta_grid,1)-1 : k_uw = k_uw

        if any((init_array[i_uw,j,k], init_array[i,j_uw,k], init_array[i,j,k_uw])) == true
            # pull value from upwind points
            val_i_uw = value_array[i_uw, j, k]
            val_j_uw = value_array[i, j_uw, k]
            val_k_uw = value_array[i, j, k_uw]

            # calculate value for given action
            val_ijk_u = finite_diff_eqn(xdot, val_i_uw, val_j_uw, val_k_uw, env)

            if val_ijk_u < val_ijk_min
                val_ijk_min = val_ijk_u
            end

            init_ijk = true
        end
    end

    if init_ijk == true
        val_ijk = val_ijk_min
    else
        val_ijk = value_array[i,j,k]
    end

    return val_ijk, init_ijk
end

function finite_diff_eqn(xdot::SVector{3, Float64}, val_i_uw::Float64, val_j_uw::Float64, val_k_uw::Float64, env::Environment)
    s1 = sign(xdot[1])
    s2 = sign(xdot[2])
    s3 = sign(xdot[3])

    num = 1.0 + s1/env.h_xy*xdot[1]*val_i_uw + s2/env.h_xy*xdot[2]*val_j_uw + s3/env.h_theta*xdot[3]*val_k_uw
    den = s1/env.h_xy*xdot[1] + s2/env.h_xy*xdot[2] + s3/env.h_theta*xdot[3]

    val_ijk = num/den

    return val_ijk
end