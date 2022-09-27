# HJB_planner_functions.jl

using GridInterpolations
using Random
using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_utils.jl")

# planner hierarchy
#   - HJB_planner()
#       - HJB_action()

# NOTE: not actually used in operation
function plan_HJB_path(x_0, dt_plan, value_array, opt_ia_array, max_steps, EoM, env, veh, sg, ag)
    if typeof(x_0) != Vector{Float64}
        x_0 = convert(Vector{Float64}, x_0)
    end

    x_k = x_0

    x_path = []  
    u_path = []
   
    step = 1
    while in_target_set(x_k, env, veh) == false && step < max_steps
        # calculate optimal action
        u_k = fast_policy(x_k, dt_plan, value_array, opt_ia_array, EoM, veh, sg, ag)  # SPEED: focus is here
        # u_k = HJB_policy(x_k, dt_plan, value_array, EoM, veh, sg, ag)

        # simulate forward one time step
        x_k1 = runge_kutta_4(x_k, u_k, dt_plan, EoM, veh, sg)    # TO-DO: needs K_sub for collision checking
        
        # store state and action at current time step
        push!(x_path, x_k)
        push!(u_path, u_k)
        
        # pass state forward to next step
        x_k = deepcopy(x_k1)
        step += 1
    end
    
    push!(x_path, x_k)

    # println("steps in HJB path: ", step)

    return x_path, u_path, step
end

# calculate one-step lookahead search at current state
function HJB_policy(x_k, dt_plan, value_array, EoM, veh, sg, ag)
    value_k1_min = Inf
    a_k_opt = ag.action_grid[1]

    for a_k in ag.action_grid
        x_k1 = runge_kutta_4(x_k, a_k, dt_plan, EoM, veh, sg)
        value_k1 = interp_value(x_k1, value_array, sg)

        if value_k1 < value_k1_min
            value_k1_min = value_k1
            a_k_opt = a_k
        end
    end

    return a_k_opt
end

# ISSUE: at certain states, returns 0 as nearest neighbor action

# TO-DO: improve action selection method
#   - goal: check as few actions as possible, while still returning a reasonable action
#   - for all possible actions, could create some ordered list of their values
#       - to be fast, want to choose action near top of list, without actually calculating list
#   - know:
#       - optimal action at each neighboring node
#       - change in value for an optimal action (equal to dt_plan)
#   - near-optimal value cutoff:
#       - in order to avoid checking whole action space for best value, need some cutoff condition to end search
#       - cutoff threshold is an estimate, could have rare case that none of the actions meet your condition
#           - in this case, will just have to search full action space and take best available
#   - ordering: 
#       - want to order indices of surrounding nodes by their distance to the state
#       - remaining action indices filled in at end, randomized
#       - no idea how to do this lol, needs to be efficient somehow

# STATUS: seems like it works...
# TO-DO: test for Dubin's car, more x_0 points

# ISSUE: struggles to get exactly into goal region
#   - think that epsilon might be too forgiving, causing it to take ugly actions
#   - when flow field converges near goal, neighboring nodes more likely to have different actions
#   - (?): does the optimal action always have Dv=-dt_plan? think interpolation/solving method get in way of this
#   - fixes:
#       - increasing epsilon (more restrictive) makes it better, but still see loops for some paths
#       - just increasing the goal works lol
#       - could tighten epsilon as you get closer to goal (use straightline dist or value knowledge)

function fast_policy(x_k, dt_plan, value_array, opt_ia_array, EoM, veh, sg, ag)
    index, weight = interpolants(sg.state_grid, x_k)
    
    ia_complete = collect(1:length(ag.action_grid))

    # @show x_k
    # @show index
    # @show weight
    ia_neighbor_srt_unq = unique(opt_ia_array[index[sortperm(weight, rev=true)]])
    deleteat!(ia_neighbor_srt_unq, ia_neighbor_srt_unq .== 0)
    
    ia_leftover_shuf_unq = shuffle(setdiff(ia_complete, opt_ia_array[index]))
    
    ia_fast_list = vcat(ia_neighbor_srt_unq, ia_leftover_shuf_unq)

    value_k = interp_value(x_k, value_array, sg)
    epsilon = 0.75 * dt_plan

    value_k1_min = Inf
    ia_min = ag.action_grid[1]

    for ia in ia_fast_list
        x_k1 = runge_kutta_4(x_k, ag.action_grid[ia], dt_plan, EoM, veh, sg)
        value_k1 = interp_value(x_k1, value_array, sg)

        Dv = value_k1 - value_k
        # println(ia)
        # println("s->a->sp (iter): ", x_k, " -> ", ag.action_grid[ia], " -> ", x_k1)
        # println("v->vp: ", value_k, " -> ", value_k1, ", Dv = ", Dv)

        # checks if tested action passes near-optimal threshold
        if value_k1 < (value_k - epsilon)
            # println("taking threshold action: ", ia, "\n")
            return ag.action_grid[ia]
        end

        # otherwise, stores best action found so far
        if value_k1 < value_k1_min
            value_k1_min = value_k1
            ia_min = ia
        end
    end

    # println("taking fully minimized action: ", ia_min, "\n")
    return ag.action_grid[ia_min]
end


# # pulling grid nodes that neighbor x_k
    # nbr_index, nbr_weight = interpolants(sg.state_grid, x_k)

    # # finding nearest neighbor action
    # nbr_index_nn = nbr_index[findmax(nbr_weight)[2]]
    # ia_pick = opt_ia_array[nbr_index_nn]

    # evaluating optimality
    # value_k = interp_value(x_k, value_array, sg)
    # value_k1 = Inf

    # if ia_pick != 0
    #     a_pick = ag.action_grid[ia_pick]
    #     x_k1 = runge_kutta_4(x_k, ag.action_grid[ia_pick], dt_plan, EoM, veh, sg)
    #     value_k1 = interp_value(x_k1, value_array, sg)
    # else
    #     a_pick = [Inf, Inf]
    #     value_k1 = Inf
    # end