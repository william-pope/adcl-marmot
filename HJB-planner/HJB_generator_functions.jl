# HJB_generator_functions.jl

using StaticArrays
using LinearAlgebra
using Rotations
using BenchmarkTools
using ProfileView

include("HJB_utils.jl")
include("dynamics_models.jl")
include("HJB_plotting.jl")

algs_path_mac = "/value_arraysers/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"

# main function to iteratively calculate continuous value function for HJB using finite difference method (STLC models only)
function solve_HJB_PDE(actions, dt_gen, dval_tol, max_steps, env::Environment, veh::Vehicle, EoM::Function, plot_growth)
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
    dval_max = Inf
    step = 1
    gs = 1

    while dval_max > dval_tol && step < max_steps
        sweep = gs_sweeps[gs]
        
        for i in sweep[1]
            for j in sweep[2]
                for k in sweep[3]
                    if target_array[i,j,k] == false && obstacle_array[i,j,k] == false
                        # # finite difference method
                        # value_array[i,j,k], init_array[i,j,k] = update_node_value_FDM(value_array, init_array, i, j, k, actions, env, veh, EoM)
                        
                        # semi-Lagrangian method
                        value_array[i,j,k] = update_node_value_SL(value_array, dt_gen, i, j, k, actions, env, veh, EoM)
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

        # plot ---
        if plot_growth == true
            theta_plot = 1/2*pi
            plot_HJB_growth(value_array, step, theta_plot, env, veh)
        end
    end

    return value_array, target_array, obstacle_array
end

# use modules
# Julia workflows

# @code_warntype

# TO-DO: may need to do collision checking along propagation path to make sure branch doesn't jump over an obstacle

# TO-DO: need more efficient way to collision check (including path) on specific actions
#   - maybe do it during initialization? will be static
#   - in basic value interpolation, one of the neighboring nodes should be +1000.0, so action won't be taken
#       - not guaranteed to always be true, checking state is more robust
#       - want to get away from using obstacle value, in which case this is probably even more important
#   - what about actions that lead to RIC? (getting into shielding stuff)
#       - is there a more systematic way of doing this than just checking a big list?

# NOTE: may need to make changes to init_array process
function update_node_value_SL(value_array, dt_gen, i::Int, j::Int, k::Int, actions, env::Environment, veh::Vehicle, EoM::Function)
    rho = veh.axle_l/tan(veh.u_phi_max)
    
    # ISSUE: value array behaves strangely when dt is small (~<= 0.1) (seems unexpected)

    x_ijk = [env.x_grid[i], env.y_grid[j], env.theta_grid[k]]

    val_p_min = Inf
    for u in actions
        x_p = runge_kutta_4(x_ijk, u, dt_gen, EoM, veh)
        val_p = interp_value(x_p, value_array, env)

        
        # if in_workspace(x_p, env, veh) == true && in_obstacle_set(x_p, env, veh) == false
        #     val_p = interp_value(x_p, value_array, env)
        # else
        #     val_p = 1000.0
        # end

        if val_p < val_p_min
            val_p_min = val_p
        end
    end

    val_ijk = dt_gen + val_p_min

    return val_ijk
end

# TO-DO: ignoring obstacle issues for now, will need to address
#   - may need to keep track of RIC (kinda like H-J reachability stuff...)
function interp_value(x::Vector{Float64}, value_array::Array{Float64, 3}, env::Environment)
    x_itp = deepcopy(x)

    # ISSUE: interpolating near boundary has same problem as near obstacles
    #   - obstacle value should not be treated as valid (but it is...)

    # adjust x,y within bounds
    x_itp[1] <= env.x_grid[1] ? x_itp[1] = env.x_grid[1]+1/16*env.h_xy : x_itp[1] = x_itp[1]
    x_itp[1] >= env.x_grid[end] ? x_itp[1] = env.x_grid[end]-1/16*env.h_xy : x_itp[1] = x_itp[1]

    x_itp[2] <= env.y_grid[1] ? x_itp[2] = env.y_grid[1]+1/16*env.h_xy : x_itp[2] = x_itp[2]
    x_itp[2] >= env.y_grid[end] ? x_itp[2] = env.y_grid[end]-1/16*env.h_xy : x_itp[2] = x_itp[2]

    # get box node indices
    i_0 = find_idx(x_itp[1], env.x_grid)
    j_0 = find_idx(x_itp[2], env.y_grid)
    k_0 = find_idx(x_itp[3], env.theta_grid)
    
    i_1 = i_0 + 1
    j_1 = j_0 + 1
    k_1 = k_0 + 1

    if k_0 == 0
        k_0 = size(env.theta_grid,1) - 1
        k_1 = 1
    end

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

                x_ijk = [x_i, y_j, theta_k]

                if in_target_set(x_ijk, env, veh) == true
                    value_array[i,j,k] = 0.0
                    init_array[i,j,k] = true
                    target_array[i,j,k] = true
                elseif in_workspace(x_ijk, env, veh) == false || in_obstacle_set(x_ijk, env, veh) == true
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