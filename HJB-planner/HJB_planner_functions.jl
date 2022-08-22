# HJB_planner_functions.jl

using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_action_functions.jl")

# planner hierarchy
#   - HJB_planner()
#       - HJB_action()
#           - value_interp()

# NOTE: not actually used in operation
function plan_HJB_path(x_0, actions, dt, value_array, obstacle_array, max_steps, EoM::Function, env::Environment, veh::Vehicle)
    if typeof(x_0) != Vector{Float64}
        x_0 = convert(Vector{Float64}, x_0)
    end

    x_k = x_0
    x_path = [x_k]

    u_path = []

    step = 0
    while in_target_set(x_k, env, veh) == false && step < max_steps
        step += 1

        # calculate optimal action
        u_k = get_HJB_action(x_k, actions, dt, value_array, obstacle_array, EoM, env, veh)  # SPEED: focus is here
        push!(u_path, u_k)

        # simulate forward one time step
        x_k1 = runge_kutta_4(x_k, u_k, dt, EoM, veh)    # needs K_sub
        push!(x_path, x_k1)

        x_k = deepcopy(x_k1)
    end

    # println("steps in HJB path: ", step)

    return x_path, u_path, step
end
