# HJB_planner_functions.jl

using GridInterpolations
using Random

include("HJB_utils.jl")

function plan_rollout_path(x_0, Dt, value_array, a_ind_opt_array, max_steps, env, veh, sg)
    x_path = []  
    x_subpath = []  
    a_path = []

    x_k = x_0
    push!(x_path, x_k)
    push!(x_subpath, x_k)
   
    step = 1
    while in_target_set(x_k, env, veh) == false && step < max_steps
        # calculate rollout action for given Dv_RC
        Dv_RC = rand([-0.5, 0.0, 0.5])
        a_k = rollout_policy(x_k, Dv_RC, Dt, value_array, veh, sg)

        # simulate forward one time step
        x_k1, x_k1_subpath = common_prop_HJB(x_k, a_k, Dt, 4)

        # store state and action at current time step
        push!(x_path, x_k1)
        for x_kk in x_k1_subpath
            push!(x_subpath, x_kk)
        end
        push!(a_path, a_k)
        
        # pass state forward to next step
        x_k = deepcopy(x_k1)
        step += 1
    end

    return x_path, x_subpath, a_path, step
end

function plan_HJB_path(x_0, Dt, value_array, a_ind_opt_array, max_steps, env, veh, sg)
    x_path = []  
    x_subpath = []  
    a_path = []

    x_k = x_0
    push!(x_path, x_k)
    push!(x_subpath, x_k)
   
    step = 1
    while in_target_set(x_k, env, veh) == false && step < max_steps
        # calculate HJB optimal action
        a_k = HJB_policy(x_k, Dt, value_array, veh, sg)

        # simulate forward one time step
        x_k1, x_k1_subpath = common_prop_HJB(x_k, a_k, Dt, 4)

        # store state and action at current time step
        push!(x_path, x_k1)
        for x_kk in x_k1_subpath
            push!(x_subpath, x_kk)
        end
        push!(a_path, a_k)
        
        # pass state forward to next step
        x_k = deepcopy(x_k1)
        step += 1
    end

    return x_path, x_subpath, a_path, step
end

# TO-DO: need to include reward in Q(s,a) optimization
# calculate one-step lookahead search at current state
function HJB_policy(x_k, Dt, value_array, veh, sg)
    val_k1_min = Inf
    a_ind_opt = 1

    actions = get_ro_action_set(x_k)

    for a_ind in eachindex(actions)
        a_k = actions[a_ind]

        x_k1, _ = common_prop_HJB(x_k, a_k, Dt, 4)
        val_k1 = interp_state_value(x_k1, value_array, sg)

        if val_k1 < val_k1_min
            val_k1_min = val_k1
            a_ind_opt = a_ind
        end
    end

    a_k_opt = actions[a_ind_opt]
    return a_k_opt
end

# use stored action grid to efficiently find near-optimal action at current state
function fast_policy(x_k, Dt, value_array, a_ind_opt_array, veh, sg)
    # gets actions from neighboring nodes
    nbr_indices, nbr_weights = interpolants(sg.state_grid, x_k)   # (!): 4 alloc (no fix)

    coord_srt = sortperm(nbr_weights, rev=true)             # (!): 2 alloc 
    
    nbr_indices_srt = view(nbr_indices, coord_srt)          # no alloc

    a_ind_neighbor_srt_unq = opt_ia_array[nbr_indices]     # (!): 1 alloc 

    unique!(ia_neighbor_srt_unq)                            # no alloc

    # assesses optimality
    value_k = interp_state_value(x_k, value_array, sg)      # (!): 5 alloc (no fix)
    epsilon = 0.75 * Dt

    actions = get_action_set(x_k)

    a_ind_opt = 1
    value_k1_min = Inf
    
    # checks neighbors first
    for a_ind in ia_neighbor_srt_unq
        if a_ind == 0
            continue
        end

        # simulates action one step forward
        x_k1, _ = common_prop_HJB(x_k, actions[a_ind], Dt, 4)
        value_k1 = interp_state_value(x_k1, value_array, sg)

        # checks if tested action passes near-optimal threshold
        if value_k1 < (value_k - epsilon)
            return actions[a_ind]
        end

        # otherwise, stores best action found so far
        if value_k1 < value_k1_min
            value_k1_min = value_k1
            a_ind_opt = a_ind
        end
    end

    # if optimal threshold hasn't been met (return), continues to check rest of actions
    a_ind_complete = collect(1:length(actions))
    a_ind_leftover_shuf_unq = shuffle(setdiff(a_ind_complete, a_ind_opt_array[nbr_indices]))

    for ia in ia_leftover_shuf_unq
        if ia == 0
            continue
        end
        
        # simulates action one step forward
        x_k1, _ = common_prop_HJB(x_k, actions[ia], Dt, 4)
        value_k1 = interp_state_value(x_k1, value_array, sg)

        # checks if tested action passes near-optimal threshold
        if value_k1 < (value_k - epsilon)
            return actions[ia]
        end

        # otherwise, stores best action found so far
        if value_k1 < value_k1_min
            value_k1_min = value_k1
            ia_min = ia
        end
    end
    
    return actions[ia_min]
end

# test by feeding random Dv and see how it does

# TO-DO: see if stored nearest-neighbor actions can be utilized
# TO-DO: write general action set optimization function for solver/planner
function rollout_policy(x_k, Dv_RC, Dt, value_array, veh, sg)
    println("\nrollout planner ---")
    println("x_k = ", x_k)
    println("Dv_RC = ", Dv_RC)

    action_HJB = HJB_policy(x_k, Dt, value_array, veh, sg)
    println("HJB action = ", action_HJB)

    # get actions from current state
    ro_actions = get_ro_action_set(x_k)

    # get value at current state (not necessary)
    val_k = interp_state_value(x_k, value_array, sg)
    println("val_k = ", val_k)

    # 1) find best phi for Dv given by reactive controller ---
    a_ind_array_RC = findall(Dv -> Dv == Dv_RC, getindex.(ro_actions, 2))

    val_best_RC = Inf 
    a_ind_best_RC = 1

    # iterate through RC actions
    for a_ind in a_ind_array_RC
        a_k = ro_actions[a_ind]

        # println("\na_ind = ", a_ind)
        # println("a_k = ", a_k)

        # propagate state and check value
        x_k1, _ = common_prop_HJB(x_k, a_k, Dt, 4)
        val_k1 = interp_state_value(x_k1, value_array, sg)

        # println("x_k1 = ", x_k1)
        # println("val_k1 = ", val_k1)

        if val_k1 < val_best_RC
            val_best_RC = val_k1
            a_ind_best_RC = a_ind
        end
    end

    println("val_best_RC = ", val_best_RC)
    println("a_ind_best_RC = ", a_ind_best_RC)
    println("action_best_RC = ", ro_actions[a_ind_best_RC])
       
    # 2) check if [Dv_RC, phi_best_RC] is a valid action in static environment ---
    infty_set_lim = 50.0
    if val_best_RC <= infty_set_lim
        a_ro = ro_actions[a_ind_best_RC]

        println("best RC action is valid")
        return a_ro
    end

    println("best RC action is not valid, moving to pure HJB")

    # 3) if RC requested Dv is not valid, then find pure HJB best action ---
    a_ind_array_no_RC = findall(Dv -> Dv != Dv_RC, getindex.(ro_actions, 2))

    val_best_HJB = Inf
    a_ind_best_HJB = 1

    # iterate through non-RC actions
    for a_ind in a_ind_array_no_RC
        a_k = ro_actions[a_ind]

        # propagate state and check value
        x_k1, _ = common_prop_HJB(x_k, a_k, Dt, 4)
        val_k1 = interp_state_value(x_k1, value_array, sg)

        if val_k1 < val_best_HJB
            val_best_HJB = val_k1
            a_ind_best_HJB = a_ind
        end
    end

    println("val_best_HJB = ", val_best_HJB)
    println("a_ind_best_HJB = ", a_ind_best_HJB)

    a_ro = ro_actions[a_ind_best_HJB]

    return a_ro
end