# HJB_planner_functions.jl

using GridInterpolations
using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_utils.jl")

# planner hierarchy
#   - HJB_planner()
#       - HJB_action()

# NOTE: not actually used in operation
function plan_HJB_path(x_0, dt_plan, value_array, action_ind_array, max_steps, EoM, env, veh, sg, ag)
    if typeof(x_0) != Vector{Float64}
        x_0 = convert(Vector{Float64}, x_0)
    end

    x_k = x_0

    x_path = []  
    u_path = []
   
    step = 1
    while in_target_set(x_k, env, veh) == false && step < max_steps
        # calculate optimal action
        u_k = fast_policy(x_k, dt_plan, value_array, action_ind_array, EoM, veh, sg, ag)  # SPEED: focus is here
       
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


# TO-DO: replace with GridInterp syntax

# calculate one-step tree search at current state
function HJB_policy(x_k, dt, value_array, EoM, veh, sg, ag)
    value_k1_min = Inf
    a_k_opt = ag.action_grid[1]

    for a_k in ag.action_grid
        x_k1 = runge_kutta_4(x_k, a_k, dt, EoM, veh, sg)
        value_k1 = interp_value(x_k1, value_array, sg)

        if value_k1 < value_k1_min
            value_k1_min = value_k1
            a_k_opt = a_k
        end
    end

    return a_k_opt
end

# ISSUE: at certain states, returns 0 for action indexing
#   - at states near obstacles/RICs, some neighboring nodes within obstacle/RIC will be invalid
#   - states in obstacles are initialized to ia = 0, so sometimes rounded interpolation retrurns 0
#   - actually a larger set of problem states, because don't want any action input
#   - (!): may be addressed by grid overapproximation that needs to happen anyway
#       - should result in agent never stepping inside a grid box that has an unsafe node in corners

# TO-DO: shouldn't be interpolating action index, since number has no meaning
#   - should be taking most common value (mode) or nearest neighbor
#   - need to include value check to show that action is near-optimal

function fast_policy(x_k, dt_plan, value_array, action_ind_array, EoM, veh, sg, ag)
    ia_opt = round(Int64, interpolate(sg.state_grid, action_ind_array, x_k))
    
    # println(x_k, " -> ", ia_opt)

    a_k_opt = ag.action_grid[ia_opt]

    return a_k_opt
end