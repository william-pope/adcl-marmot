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
                        value_array[i,j,k], init_array[i,j,k] = update_node_value(value_array, init_array, i, j, k, actions, env, veh, EoM)
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

            p_k = heatmap(env.x_grid, env.y_grid, transpose(value_array[:,:,k_plot]), clim=(0,10),
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
            plot_polygon(p_k, env.obstacle_array_vec[1], 3, :red, "Obstacle")
            for obstacle_array in env.obstacle_array_vec
                plot_polygon(p_k, obstacle_array, 3, :red, "")
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

function update_node_value(value_array, init_array, i::Int, j::Int, k::Int, actions, env::Environment, veh::Vehicle, EoM::Function)   
    x1 = env.x_grid[i]
    x2 = env.y_grid[j]
    x3 = env.theta_grid[k]

    x = SVector{3, Float64}(x1, x2, x3)

    # compute val_ijk update
    value_array_min = Inf
    init_ijk = false

    # code in this loop gets iterated the most
    for ia in eachindex(actions)
        xdot = EoM(x, actions[ia], veh)
    
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
            val_a_ijk = finite_diff_eqn(xdot, val_i_uw, val_j_uw, val_k_uw, env)

            if val_a_ijk < value_array_min
                value_array_min = val_a_ijk
            end

            init_ijk = true
        end
    end

    if init_ijk == true
        value_array_ijk = value_array_min
    else
        value_array_ijk = value_array[i,j,k]
    end

    return value_array_ijk, init_ijk
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