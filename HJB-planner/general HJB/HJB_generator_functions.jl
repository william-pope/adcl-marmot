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
function solve_HJB_PDE(env, veh, EoM, sg, actions, dt_solve, dval_tol, max_solve_steps, plot_growth_flag)
    # initialize value_array, value_array, init_array
    value_array, init_array, target_array, obstacle_array = initialize_value_array(sg, env, veh)
    value_array_old = deepcopy(value_array)

    return value_array

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
                        value_array[i,j,k] = update_node_value_SL(value_array, dt_solve, i, j, k, actions, env, veh, EoM)
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

    return value_array
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
function update_node_value_SL(value_array, dt_solve, i::Int, j::Int, k::Int, actions, env::Environment, veh::Vehicle, EoM::Function)    
    # ISSUE: value array behaves strangely when dt is small (~<= 0.1) (seems unexpected)

    x_ijk = [env.x_grid[i], env.y_grid[j], env.theta_grid[k]]

    qval_min = Inf
    for u in actions
        cost_p = get_cost(x_ijk, u, dt_solve)

        x_p = runge_kutta_4(x_ijk, u, dt_solve, EoM, veh)
        val_p = interp_value(x_p, value_array, env)

        qval_u = cost_p + val_p

        # # implements obstacle checking for exact propagated state
        # if in_workspace(x_p, env, veh) == true && in_obstacle_set(x_p, env, veh) == false
        #     val_p = interp_value(x_p, value_array, env)
        # else
        #     val_p = 1000.0
        # end

        if qval_u < qval_min
            qval_min = qval_u
        end
    end

    val_ijk = qval_min

    return val_ijk
end

function get_cost(x_k, u_k, dt_solve)
    cost_k = dt_solve

    if x_k[1] >= 4.0
        cost_k = 1/2*cost_k
    end

    return cost_k
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
function initialize_value_array(sg, env, veh)
    value_array = zeros(Float64, sg.grid_size_array...)
    init_array = zeros(Bool, sg.grid_size_array...)
    target_array = zeros(Bool, sg.grid_size_array...)
    obstacle_array = zeros(Bool, sg.grid_size_array...)

    # TO-DO: iterating over rows is slowe, need to replace with more performant method (especially in solver)
    for i in 1:size(sg.grid_idx_matrix,1)
        grid_idx = sg.grid_idx_matrix[i,:]
        x_ijk = idx_to_state(grid_idx, sg)

        # TO-DO: need to update set checkers to use DomainSets.jl (or revert to old method)
        if in_workspace(x_ijk, env, veh) == false || in_obstacle_set(x_ijk, env, veh) == true
            value_array[grid_idx...] = 1000.0
            init_array[grid_idx...] = false
            obstacle_array[grid_idx...] = true
        elseif in_target_set(x_ijk, env, veh) == true
            value_array[grid_idx...] = 0.0
            init_array[grid_idx...] = true
            target_array[grid_idx...] = true
        else
            value_array[grid_idx...] = 100.0
            init_array[grid_idx...] = false
            target_array[grid_idx...] = false
            obstacle_array[grid_idx...] = false
        end
    end

    return value_array, init_array, target_array, obstacle_array
end

function idx_to_state(grid_idx, sg)
    x_ijk = zeros(Float64, sg.state_dim)

    for (d, i) in enumerate(grid_idx)
        x_ijk[d] = sg.grid_array[d][i] 
    end

    return x_ijk
end