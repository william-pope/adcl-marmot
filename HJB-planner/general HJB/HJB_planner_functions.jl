# HJB_planner_functions.jl

using GridInterpolations
using Random

include("HJB_utils.jl")

# NOTE: not actually used in operation
function plan_HJB_path(x_0, Dt, value_array, opt_ia_array, max_steps, EoM, env, veh, sg)
    x_path = []  
    x_subpath = []  
    a_path = []

    x_k = x_0
    push!(x_path, x_k)
    push!(x_subpath, x_k)
   
    step = 1
    while in_target_set(x_k, env, veh) == false && step < max_steps
        # calculate optimal action
        # a_k = fast_policy(x_k, Dt, value_array, opt_ia_array, EoM, veh, sg)
        a_k = HJB_policy(x_k, Dt, value_array, EoM, veh, sg)

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

# calculate one-step lookahead search at current state
function HJB_policy(x_k, Dt, value_array, EoM, veh, sg)
    val_k1_min = Inf
    ia_k_opt = 1

    actions = get_action_set(x_k)

    for ia_k in eachindex(actions)
        a_k = actions[ia_k]

        x_k1, _ = common_prop_HJB(x_k, a_k, Dt, 4)
        val_k1 = interp_value(x_k1, value_array, sg)

        if val_k1 < val_k1_min
            val_k1_min = val_k1
            ia_k_opt = ia_k
        end
    end

    a_k_opt = actions[ia_k_opt]
    return a_k_opt
end

# use stored action grid to efficiently find near-optimal action at current state
function fast_policy(x_k, Dt, value_array, opt_ia_array, EoM, veh, sg)
    # gets actions from neighboring nodes
    nbr_indices, nbr_weights = interpolants(sg.state_grid, x_k)   # (!): 4 alloc (no fix)

    coord_srt = sortperm(nbr_weights, rev=true)             # (!): 2 alloc 
    
    nbr_indices_srt = view(nbr_indices, coord_srt)          # no alloc

    ia_neighbor_srt_unq = opt_ia_array[nbr_indices]     # (!): 1 alloc 

    unique!(ia_neighbor_srt_unq)                            # no alloc

    # assesses optimality
    value_k = interp_value(x_k, value_array, sg)      # (!): 5 alloc (no fix)
    epsilon = 0.75 * Dt

    actions = get_action_set(x_k)

    ia_min = 1
    value_k1_min = Inf
    
    # checks neighbors first
    for ia in ia_neighbor_srt_unq
        if ia == 0
            continue
        end

        # simulates action one step forward
        x_k1, _ = common_prop_HJB(x_k, actions[ia], Dt, 4)
        value_k1 = interp_value(x_k1, value_array, sg)

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

    # if optimal threshold hasn't been met (return), continues to check rest of actions
    ia_complete = collect(1:length(actions))
    ia_leftover_shuf_unq = shuffle(setdiff(ia_complete, opt_ia_array[nbr_indices]))

    for ia in ia_leftover_shuf_unq
        if ia == 0
            continue
        end
        
        # simulates action one step forward
        x_k1, _ = common_prop_HJB(x_k, actions[ia], Dt, 4)
        value_k1 = interp_value(x_k1, value_array, sg)

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