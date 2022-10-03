# HJB_action_functions.jl

# ISSUE: one-step planner might have same issues with value interpolation
#   - possible that all 'forward' actions pick up a grid node within the obstacle
#   - due to small time steps, all x_k1 points probably in same grid cell

#   - also they're all returning val_k1=0 lol...
#       - this was added to deal with obstacles, obviously need smarter approach
#       - need to consider obs_k1
#           - still going to break at certain states (when all obs/when fwd obs)

# new idea: add a small boundary around each obstacle, use same value function
# - vehicle cannot choose actions that cross the boundary
# - vehicle can use value points behind boundary to inform gradient/tree for other actions
# - objective: can never get close enough to obstacle to break value interpolation

# new idea: use value extrapolation from nearby grids that are obstacle free

# calculate one-step tree search at current state
function get_HJB_action(x_k, actions, dt, value_array, obstacle_array, EoM::Function, env::Environment, veh::Vehicle)
    val_k1_min = Inf
    u_k_opt = actions[1]
    for u_k in actions
        x_k1 = runge_kutta_4(x_k, u_k, dt, EoM, veh)
        val_k1 = interp_value(x_k1, value_array, env)

        # println("u_k: ", u_k, ", x_k1: ", x_k1, ", val_k1: ", val_k1)

        if val_k1 < val_k1_min
            val_k1_min = val_k1
            u_k_opt = u_k
        end
    end

    # println("x_k: ", x_k)
    # println("u_opt: ", u_opt)
    return u_k_opt
end