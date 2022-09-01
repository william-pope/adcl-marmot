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

# STATUS: GridInterp.jl looks good to go, just need to rework around it
#   - still need to figure out specifics for Gauss-Seidel flipping

# main function to iteratively calculate continuous value function for HJB using finite difference method (STLC models only)
function solve_HJB_PDE(env, veh, EoM, sg, actions, dt_solve, dval_tol, max_solve_steps, plot_growth_flag)
    # initialize value_array, value_array, init_array
    value_array, init_array, target_array, obstacle_array = initialize_value_array(sg, env, veh)
    value_array_old = deepcopy(value_array)

    return value_array

    # TO-DO: rewrite Gauss-Seidel mechanism to use new matrix
    # - should be similar to initialize() use, but need to flip columns in the matrix to do different sweeps

    # main function loop
    dval_max = Inf
    solve_step = 1
    gs_step = 1

    while dval_max > dval_tol && solve_step < max_solve_steps
        # flip rows in idx matrix according to Gauss-Seidel scheme
        gs_idx_matrix = zeros(Float64, size(sg.grid_idx_matrix))
        for d in sg.gs_sweep_matrix[:,gs_step]
            if d == 1
                gs_idx_matrix[d,:] = reverse(sg.gs_sweep_matrix[d,:])
            else
                gs_idx_matrix[d,:] = sg.gs_sweep_matrix[d,:]
            end
        end

        for j in 1:sg.num_nodes
            grid_idx = sg.gs_idx_matrix[:,j]
            x_ijk = idx_to_state(grid_idx, sg)

            if target_array[ind_s] == false && obstacle_array[ind_s] == false
                # # finite difference method
                # value_array[i,j,k], init_array[i,j,k] = update_node_value_FDM(value_array, init_array, i, j, k, actions, env, veh, EoM)
                
                # semi-Lagrangian method
                value_array[i,j,k] = update_node_value_SL(x_ijk, value_array, actions, dt_solve, EoM, env, veh, sg)
            end
        end
        
        # TO-DO: see if deepcopy() can be removed
        value_array[:,:,end] = deepcopy(value_array[:,:,1])
        init_array[:,:,end] = deepcopy(init_array[:,:,1])

        # compare value_array and value_array to check convergence
        dval = value_array - value_array_old
        dval_max = maximum(abs.(dval))

        # TO-DO: see if deepcopy() can be removed
        value_array_old = deepcopy(value_array)

        println("step: ", step, ", dval_max = ", dval_max)

        if gs_step == 2^sg.state_dim
            gs_step = 1
        else
            gs_step += 1
        end

        solve_step += 1

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
function update_node_value_SL(x_ijk, value_array, actions, dt_solve, EoM, env, veh, sg) 
    # ISSUE: value array behaves strangely when dt is small (~<= 0.1) (seems unexpected)

    qval_min = Inf
    for u in actions
        cost_p = get_cost(x_ijk, u, dt_solve)

        x_p = runge_kutta_4(x_ijk, u, dt_solve, EoM, veh, sg)
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

    return cost_k
end

# NOTE: may still need this function to wrap around interpolate() to do obstacle grid overestimation stuff
# TO-DO: ignoring obstacle issues for now, will need to address
#   - may need to keep track of RIC (kinda like H-J reachability stuff...) in order to track "invalid" nodes
#   - don't need to check if neighbors are invalid every time, can store during initialization if node borders an obstacle
function interp_value(x::Vector{Float64}, value_array::Array{Float64, 3}, env::Environment)
    x_itp = deepcopy(x)

    # ISSUE: interpolating near boundary has same problem as near obstacles
    #   - obstacle value should not be treated as valid (but it is...)

    # adjust x,y within bounds
    x_itp[1] <= env.x_grid[1] ? x_itp[1] = env.x_grid[1]+1/16*env.h_xy : x_itp[1] = x_itp[1]
    x_itp[1] >= env.x_grid[end] ? x_itp[1] = env.x_grid[end]-1/16*env.h_xy : x_itp[1] = x_itp[1]

    x_itp[2] <= env.y_grid[1] ? x_itp[2] = env.y_grid[1]+1/16*env.h_xy : x_itp[2] = x_itp[2]
    x_itp[2] >= env.y_grid[end] ? x_itp[2] = env.y_grid[end]-1/16*env.h_xy : x_itp[2] = x_itp[2]

    return val_itp
end

# initialize value approximations
function initialize_value_array(sg, env, veh)
    # TO-DO: need to figure out indexing for data arrays
    #   - if data array has to be 1-d, can write some equation as f(ind) to calc 1-d ind
    value_array = zeros(Float64, length(sg.state_grid))
    init_array = zeros(Bool, length(sg.state_grid))
    target_array = zeros(Bool, length(sg.state_grid))
    obstacle_array = zeros(Bool, length(sg.state_grid))

    for ind_m in sg.ind_gs_array[1]
        x = sg.state_grid[ind_m...]
        ind_s = multi2single_ind(ind_m, sg)

        # println(ind_m, " -> ", ind_s, " -> ", x)

        if in_workspace(x, env, veh) == false || in_obstacle_set(x, env, veh) == true
            value_array[ind_s] = 1000.0
            init_array[ind_s] = false
            obstacle_array[ind_s] = true
        elseif in_target_set(x, env, veh) == true
            value_array[ind_s] = 0.0
            init_array[ind_s] = true
            target_array[ind_s] = true
        else
            value_array[ind_s] = 100.0    # TO-DO: make 1000.0
            init_array[ind_s] = false
            target_array[ind_s] = false
            obstacle_array[ind_s] = false
        end
    end

    return value_array, init_array, target_array, obstacle_array
end