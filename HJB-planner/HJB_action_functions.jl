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
        val_k1, obs_k1 = interp_value_act(x_k1, value_array, obstacle_array, env)

        # println("u_k: ", u_k, ", x_k1: ", x_k1, ", val_k1: ", val_k1, ", obs_k1: ", obs_k1)

        if val_k1 < val_k1_min && obs_k1 == false
            val_k1_min = val_k1
            u_k_opt = u_k
        end
    end

    # println("x_k: ", x_k)
    # println("u_opt: ", u_opt)
    return u_k_opt
end

# may produce asymmetry due to interpolation index selection
function interp_value_act(x::Vector{Float64}, U::Array{Float64, 3}, O::Array{Bool, 3}, env::Environment)
    # TO-DO: need to address boundary and deepcopy(x)

    # x_itp = deepcopy(x)     # SPEED: big impact on allocations

    # # ?: if weirdly sharp curves still happening, check if this is reassigning state
    # # map out-of-bounds xy states to boundary
    # if x_itp[1] > env.x_grid[end-1]
    #     x_itp[1] = env.x_grid[end-1]
    # elseif x_itp[1] < env.x_grid[2]
    #     x_itp[1] = env.x_grid[2]
    # end

    # if x_itp[2] > env.y_grid[end-1]
    #     x_itp[2] = env.y_grid[end-1]
    # elseif x_itp[2] < env.y_grid[2]
    #     x_itp[2] = env.y_grid[2]
    # end
    
    x_itp = x

    # finds indices of neighboring grid nodes
    #   - made up of the 8 grid points that surround target state
    #   - need to test all 8 points for obstacle collisions
    #   - if any of 8 is in obstacle, whole point gets thrown out
    #   - can't use value from obstacle, since it isn't part of flow field

    i_0 = find_idx(x_itp[1], env.x_grid)
    j_0 = find_idx(x_itp[2], env.y_grid)
    k_0 = find_idx(x_itp[3], env.theta_grid)

    i_1 = i_0 + 1
    j_1 = j_0 + 1
    k_0 == size(env.theta_grid,1) ? k_1 = 1 : k_1 = k_0 + 1

    # have indices for interpolation grid, need to check for obstacles
    obs_itp = false

    # for i in (i_0, i_1)
    #     for j in (j_0, j_1)
    #         for k in (k_0, k_1)
    #             # println("[i,j,k]: ", [i,j,k], ", [x,y,th]: ", [env.x_grid[i], env.y_grid[j], env.theta_grid[k]], ", U[i,j,k] = ", U[i,j,k])

    #             # if any nodes used for interp are in obstacle, whole interpolation for "x_1p" is thrown out
    #             if O[i,j,k] == true
    #                 obs_itp = true
    #                 u_itp = 0.0

    #                 return u_itp, obs_itp
    #             end
    #         end
    #     end
    # end

    x_0 = env.x_grid[i_0]
    y_0 = env.y_grid[j_0]
    theta_0 = env.theta_grid[k_0]

    x_1 = env.x_grid[i_1]   
    y_1 = env.y_grid[j_1]
    theta_1 = env.theta_grid[k_1]

    x_d = (x_itp[1] - x_0)/(x_1 - x_0)
    y_d = (x_itp[2] - y_0)/(y_1 - y_0)
    theta_d = (x_itp[3] - theta_0)/(theta_1 - theta_0)

    u_000 = U[i_0, j_0, k_0]
    u_100 = U[i_1, j_0, k_0]
    u_010 = U[i_0, j_1, k_0]
    u_110 = U[i_1, j_1, k_0]
    u_001 = U[i_0, j_0, k_1]
    u_101 = U[i_1, j_0, k_1]
    u_011 = U[i_0, j_1, k_1]
    u_111 = U[i_1, j_1, k_1]

    u_00 = u_000*(1 - x_d) + u_100*x_d
    u_01 = u_001*(1 - x_d) + u_101*x_d
    u_10 = u_010*(1 - x_d) + u_110*x_d
    u_11 = u_011*(1 - x_d) + u_111*x_d

    u_0 = u_00*(1 - y_d) + u_10*y_d
    u_1 = u_01*(1 - y_d) + u_11*y_d

    u_itp = u_0*(1 - theta_d) + u_1*theta_d

    return u_itp, obs_itp
end