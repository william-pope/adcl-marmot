# HJB_solver_functions.jl

using StaticArrays

include("HJB_utils.jl")
include("HJB_plotting.jl")

# main function to iteratively calculate HJB value function
function solve_HJB_PDE(env, veh, sg, Dt, Dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)
    # initialize data arrays
    value_array, a_ind_opt_array, set_array = initialize_value_array(sg, env, veh)

    num_gs_sweeps = 2^dimensions(sg.state_grid)

    # main function loop
    Dval_max = Inf
    solve_step = 1
    gs_step = 1

    while solve_step <= max_solve_steps
        Dval_max = 0.0

        for ind_m in sg.ind_gs_array[gs_step]
            ind_s = multi2single_ind(ind_m, sg)

            # if the node is in free space, update its value
            if set_array[ind_s] == 2
                x = sg.state_list_static[ind_s]
                
                # store previous value
                v_kn1 = value_array[ind_s]
                
                # calculate new value
                value_array[ind_s], a_ind_opt_array[ind_s] = update_node_value(x, value_array, Dt, env, veh, sg)
                
                # compare old and new values, update largest change in value
                v_k = value_array[ind_s]
                Dval = abs(v_k - v_kn1)

                if Dval > Dval_max
                    Dval_max = Dval
                end
            end
        end

        println("solve_step: ", solve_step, ", gs_step: ", gs_step, ", Dval_max = ", Dval_max)

        # check if convergence condition is met
        if Dval_max <= Dval_tol
            break
        end

        # update step counters
        if gs_step == num_gs_sweeps
            gs_step = 1
        else
            gs_step += 1
        end

        solve_step += 1

        # plot ---
        if plot_growth_flag == 1
            plot_HJB_growth(value_array, heatmap_clim, solve_step, env, veh)
        end
    end

    return value_array, a_ind_opt_array, set_array
end

# TO-DO: may need to do collision-checking along propagation subpath
function update_node_value(x, value_array, Dt, env, veh, sg) 
    qval_min = Inf
    a_ind_opt = 1

    actions = get_ro_action_set(x)

    for a_ind in eachindex(actions)
        a = actions[a_ind]

        cost_p = get_cost(x, a, Dt)

        x_p, x_p_subpath = common_prop_HJB(x, a, Dt, 4)

        val_p = interp_subpath_value(x_p, x_p_subpath, value_array, sg)

        qval_a = cost_p + val_p

        if qval_a < qval_min
            qval_min = qval_a
            a_ind_opt = a_ind
        end
    end

    val_ijk = qval_min
   
    return val_ijk, a_ind_opt
end

# initialize arrays
function initialize_value_array(sg, env, veh)
    value_array = zeros(Float64, length(sg.state_grid))
    a_ind_opt_array = zeros(Int, length(sg.state_grid))
    set_array = zeros(Int, length(sg.state_grid))

    for ind_m in sg.ind_gs_array[1]
        ind_s = multi2single_ind(ind_m, sg)
        x = sg.state_list_static[ind_s]

        if in_workspace(x, env, veh) == false || in_obstacle_set(x, env, veh) == true
            value_array[ind_s] = 1e5
            set_array[ind_s] = 0
        
        elseif in_target_set(x, env, veh) == true
            value_array[ind_s] = 0.0
            set_array[ind_s] = 1
        
        else
            value_array[ind_s] = 1e5
            set_array[ind_s] = 2
        end
    end

    return value_array, a_ind_opt_array, set_array
end