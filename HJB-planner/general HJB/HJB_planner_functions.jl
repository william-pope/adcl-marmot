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
    ro_actions = get_ro_action_set(x_k)

    a_ind_array = collect(1:length(ro_actions))
    a_ind_opt, _ = optimize_action(x_k, a_ind_array, ro_actions, value_array, Dt, sg)

    a_k_opt = ro_actions[a_ind_opt]

    return a_k_opt
end

# TO-DO: see if stored nearest-neighbor actions can be utilized
function rollout_policy(x_k, Dv_RC, Dt, value_array, veh, sg)    
    # get actions for current state
    ro_actions = get_ro_action_set(x_k)

    # 1) find best phi for Dv given by reactive controller ---
    a_ind_array_RC = findall(Dv -> Dv == Dv_RC, getindex.(ro_actions, 2))
    a_ind_best_RC, val_best_RC = optimize_action(x_k, a_ind_array_RC, ro_actions, value_array, Dt, sg)

    # 2) check if [Dv_RC, phi_best_RC] is a valid action in static environment ---
    infty_set_lim = 50.0
    if val_best_RC <= infty_set_lim
        a_ro = ro_actions[a_ind_best_RC]
 
        return a_ro
    end

    # 3) if RC requested Dv is not valid, then find pure HJB best action ---
    a_ind_array_no_RC = findall(Dv -> Dv != Dv_RC, getindex.(ro_actions, 2))
    a_ind_best_HJB, val_best_HJB = optimize_action(x_k, a_ind_array_no_RC, ro_actions, value_array, Dt, sg)

    a_ro = ro_actions[a_ind_best_HJB]

    return a_ro
end

# println("\nrollout planner ---")
# println("x_k = ", x_k)
# println("Dv_RC = ", Dv_RC)

# action_HJB = HJB_policy(x_k, Dt, value_array, veh, sg)
# println("HJB action = ", action_HJB)

# # get value at current state (not necessary)
# val_k = interp_state_value(x_k, value_array, sg)
# println("val_k = ", val_k)

# println("val_best_RC = ", val_best_RC)
# println("a_ind_best_RC = ", a_ind_best_RC)
# println("action_best_RC = ", ro_actions[a_ind_best_RC])

# println("best RC action is valid")

# println("best RC action is not valid, moving to pure HJB")

# println("val_best_HJB = ", val_best_HJB)
# println("a_ind_best_HJB = ", a_ind_best_HJB)