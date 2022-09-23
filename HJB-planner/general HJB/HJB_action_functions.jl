# HJB_action_functions.jl

# TO-DO: replace with GridInterp syntax

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