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

# TO-DO: need to store optimal actions during solve

# main function to iteratively calculate HJB value function
function solve_HJB_PDE(env, veh, EoM, sg, ag, dt_solve, dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)
    # initialize value_array, init_array, ...
    value_array, action_ind_array, init_array, target_array, obstacle_array = initialize_value_array(sg, env, veh)
    value_array_old = deepcopy(value_array)

    # main function loop
    dval_max = Inf
    solve_step = 1
    gs_step = 1

    # SPEED:
    #   - with deepcopy(value_array) -> 39.7 s, 292,707,617 alloc
    #   - with old[ind_s] = new[ind_s] -> 39.7 s, 293,412,717 alloc (basically no change)
    #   - need to make data_array static

    # dval_max > dval_tol && 
    while solve_step <= max_solve_steps
        for ind_m in sg.ind_gs_array[gs_step]
            x = sg.state_grid[ind_m...]
            ind_s = multi2single_ind(ind_m, sg)

            if target_array[ind_s] == false && obstacle_array[ind_s] == false
                # value_array_old[ind_s] = value_array[ind_s]
                value_array[ind_s], action_ind_array[ind_s] = update_node_value_SL(x, value_array, dt_solve, EoM, env, veh, sg, ag)
            end
        end
        
        # TO-DO: see if this is needed for angle wrap
        # value_array[:,:,end] = value_array[:,:,1]
        # init_array[:,:,end] = deepcopy(init_array[:,:,1])

        # compare value_array and value_array to check convergence
        dval = value_array - value_array_old
        dval_max = maximum(abs.(dval))

        # TO-DO: see if deepcopy() can be replaced
        value_array_old = deepcopy(value_array)

        println("step: ", solve_step, ", dval_max = ", dval_max)

        if gs_step == 2^dimensions(sg.state_grid)
            gs_step = 1
        else
            gs_step += 1
        end

        solve_step += 1

        # plot ---
        if plot_growth_flag == true
            k_plot = 19
            plot_HJB_growth(value_array, heatmap_clim, step, k_plot, env, veh)
        end
    end

    return value_array, action_ind_array
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
function update_node_value_SL(x, value_array, dt_solve, EoM, env, veh, sg, ag) 
    # ISSUE: value array behaves strangely when dt is small (~<= 0.1) (seems unexpected)

    qval_min = Inf
    ia_opt_ijk = 1

    for ia in 1:length(ag.action_grid)
        a = ag.action_grid[ia]

        cost_p = get_cost(x, a, dt_solve)

        x_p = runge_kutta_4(x, a, dt_solve, EoM, veh, sg)
        val_p = interp_value(x_p, value_array, sg)

        qval_a = cost_p + val_p

        if qval_a < qval_min
            qval_min = qval_a
            ia_opt_ijk = ia
        end
    end

    val_ijk = qval_min

    return val_ijk, ia_opt_ijk
end

function get_cost(x_k, a_k, dt_solve)
    cost_k = dt_solve

    return cost_k
end

# NOTE: may still need this function to wrap around interpolate() to do obstacle grid overestimation stuff
# TO-DO: ignoring obstacle issues for now, will need to address
#   - may need to keep track of RIC (kinda like H-J reachability stuff...) in order to track "invalid" nodes
#   - don't need to check if neighbors are invalid every time, can store during initialization if node borders an obstacle
function interp_value(x, value_array, sg)
    # check if current state is within state space
    for d in eachindex(x)
        if x[d] < sg.state_grid.cutPoints[d][1] || x[d] > sg.state_grid.cutPoints[d][end]
            val_itp = 1000.0
            return val_itp
        end
    end

    val_itp = interpolate(sg.state_grid, value_array, x)

    return val_itp
end

# initialize value approximations
function initialize_value_array(sg, env, veh)
    # TO-DO: need to figure out indexing for data arrays
    #   - if data array has to be 1-d, can write some equation as f(ind) to calc 1-d ind
    value_array = zeros(Float64, length(sg.state_grid))
    action_ind_array = zeros(Int, length(sg.state_grid))
    init_array = zeros(Bool, length(sg.state_grid))
    target_array = zeros(Bool, length(sg.state_grid))
    obstacle_array = zeros(Bool, length(sg.state_grid))

    for ind_m in sg.ind_gs_array[1]
        x = sg.state_grid[ind_m...]
        ind_s = multi2single_ind(ind_m, sg)

        # println(ind_m, " -> ", ind_s, " -> ", x)

        # SPEED: in_obstacle_set() is horrible for performance
        if in_workspace(x, env, veh) == false || in_obstacle_set(x, env, veh) == true
            value_array[ind_s] = 1e5
            init_array[ind_s] = false
            obstacle_array[ind_s] = true
        
        elseif in_target_set(x, env, veh) == true
            value_array[ind_s] = 0.0
            init_array[ind_s] = true
            target_array[ind_s] = true
        
        else
            value_array[ind_s] = 100.0   # TO-DO: make 1000.0
            init_array[ind_s] = false
            target_array[ind_s] = false
            obstacle_array[ind_s] = false
        end
    end

    return value_array, action_ind_array, init_array, target_array, obstacle_array
end