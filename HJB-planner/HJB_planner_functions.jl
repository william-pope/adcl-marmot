

using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")

# planner hierarchy
#   - HJB_planner()
#       - HJB_action()
#           - value_gradient()
#               - value_interp()

# planner hierarchy (alternate)
#   - HJB_planner()
#       - HJB_action()
#           - value_interp()

# NOTE: not actually used in operation
function HJB_planner(x_0, U, dt, plan_steps, A, O, EoM::Function, env::Environment, veh::Vehicle)
    if typeof(x_0) != Vector{Float64}
        x_0 = convert(Vector{Float64}, x_0)
    end

    x_k = x_0
    x_path = [x_k]  

    u_path = []

    step = 0
    while in_target_set(x_k, env, veh) == false && step < plan_steps
        step += 1

        # calculate optimal action
        u_k = HJB_action(x_k, U, A, O, dt, EoM, env, veh)  # SPEED: focus is here
        push!(u_path, u_k)

        # simulate forward one time step
        x_k1 = runge_kutta_4(x_k, u_k, dt, EoM, veh)    # needs K_sub
        push!(x_path, x_k1)

        x_k = deepcopy(x_k1)
    end

    # println("steps in HJB path: ", step)

    return x_path, u_path, step
end

# calculates optimal action as a function of state
function HJB_action(x_k, U, A, O, dt, EoM::Function, env::Environment, veh::Vehicle)
    # TO-DO: make sure this isn't causing any issues, conflicts with other boundary mechanisms
    # adjust position when near edge for boundary issues
    if abs(x_k[1] - env.x_grid[1]) < 2*env.h_xy
        x_k[1] = env.x_grid[3]
    elseif abs(x_k[1] - env.x_grid[end]) < 2*env.h_xy
        x_k[1] = env.x_grid[end-2]
    end

    if abs(x_k[2] - env.y_grid[1]) < 2*env.h_xy
        x_k[2] = env.y_grid[3]
    elseif abs(x_k[2] - env.y_grid[end]) < 2*env.h_xy
        x_k[2] = env.y_grid[end-2]
    end

    x_k[3] > pi ? x_k[3] = x_k[3] - 2*pi : x_k[3] = x_k[3]

    println("\n--- ---\nx_k: ", x_k)
    
    # # calculate gradient at current state
    # dV = value_gradient(x_k, U, O, env)       # SPEED: big focus

    # vdot_min = Inf
    # a_opt = A[1]
    # for i in 1:size(A,1)
    #     xs = SVector{3, Float64}(x_k)     # TO-DO: clean up variable types
    #     xdot = car_EoM(xs, A[i], veh)

    #     vdot = dot(dV, xdot)

    #     if vdot < vdot_min
    #         vdot_min = vdot
    #         a_opt = A[i]
    #     end
    # end

    # ISSUE: one-step planner might have same issues with value interpolation
    #   - possible that all 'forward' actions pick up a grid node within the obstacle
    #   - due to small time steps, all x_k1 points probably in same grid cell

    #   - also they're all returning v_k1=0 lol...
    #       - this was added to deal with obstacles, obviously need smarter approach
    #       - need to consider obs_k1
    #           - still going to break at certain states (when all obs/when fwd obs)

    # new idea: add a small boundary around each obstacle, use same value function
    # - vehicle cannot choose actions that cross the boundary
    # - vehicle can use value points behind boundary to inform gradient/tree for other actions
    # - objective: can never get close enough to obstacle to break value interpolation

    # new idea: use value extrapolation from nearby grids that are obstacle free

    # calculate one-step tree search at current state
    v_k1_min = Inf
    a_opt = A[1]
    for i in 1:size(A,1)
        x_k1 = runge_kutta_4(x_k, A[i], dt, EoM, veh)
        v_k1, obs_k1 = value_interp(x_k1, U, O, env)

        println("a_k: ", A[i], ", x_k1: ", x_k1, ", v_k1: ", v_k1, ", obs_k1: ", obs_k1)

        if v_k1 < v_k1_min && obs_k1 == false
            v_k1_min = v_k1
            a_opt = A[i]
        end
    end

    # println("x_k: ", x_k)
    println("a_opt: ", a_opt)
    return a_opt
end

# ISSUE: dealing with obstacles during gradient approximation
#   - unsure what to do when one of interpolation nodes is in an obstacle
#   - depending on grid size, interp node could be somewhat far away from actual point of interest
#   - need to somehow calculate a gradient without using affected nodes
#   - gradient process is only a function of the state, not related to actions being tested
#   - obstacles are not a part of the flow field, treat like blank regions
#   - partials can be set to 0 if needed, just implies a flat gradient along that axis
#       - probably more accurate to switch to single point difference method, but will be a little harder

function value_gradient(x::Vector{Float64}, U::Array{Float64, 3}, O::Array{Bool, 3}, env::Environment)
    # println(x)

    # ISSUE: chattering still happening near obstacles
    #   - in this case, think neighbor in x and theta are obstacle states
    #   - actually negative y is too (doesn't look like it should...)
    #   - obstacles on all three axes -> frozen
    #   - single point gradient will fix this

    #   - each axis gradient is independent, can look at each in a for loop
    #   - logic:
    #       - if both free, use center difference method with both points
    #       - if one point in obs, use single diff method with other point and x
    #       - if both in obs, set du_dy = 0.0

    #   - think issue is due to both outlier point and central point touching obstacle
    
    # interpolate value at surrounding points
    v_x, obs_x = value_interp(x, U, O, env)
    v_1p, obs_1p = value_interp(x + [env.h_xy, 0.0, 0.0], U, O, env)
    v_1n, obs_1n = value_interp(x - [env.h_xy, 0.0, 0.0], U, O, env)
    v_2p, obs_2p = value_interp(x + [0.0, env.h_xy, 0.0], U, O, env)
    v_2n, obs_2n = value_interp(x - [0.0, env.h_xy, 0.0], U, O, env)
    v_3p, obs_3p = value_interp(x + [0.0, 0.0, env.h_theta], U, O, env)
    v_3n, obs_3n = value_interp(x - [0.0, 0.0, env.h_theta], U, O, env)

    # ISSUE: interpolation at center point can be in obstacle
    #   - occurs just after vehicle rounds corner of circle (9 o'clock position)
    #   - y_neg in obstacle seems strange, since behind vehicle
    #   - this might just be what happens on a diagonal face
    #   - happens when making sharp turns around obstacles (need more theta steps?)
    #   - is step size just too large? seems like we're already at the limit

    # calculating du_dx
    if obs_1p == false && obs_1n == false
        du_dx = (v_1p - v_1n)/(2*env.h_xy)
    elseif obs_1p == true && obs_1n == false
        du_dx = (v_x - v_1n)/env.h_xy
    elseif obs_1p == false && obs_1n == true
        du_dx = (v_1p - v_x)/env.h_xy
    else
        du_dx = 0.0
    end

    # calculating du_dy
    if obs_2p == false && obs_2n == false
        du_dy = (v_2p - v_2n)/(2*env.h_xy)
    elseif obs_2p == true && obs_2n == false
        du_dy = (v_x - v_2n)/env.h_xy
    elseif obs_2p == false && obs_2n == true
        du_dy = (v_2p - v_x)/env.h_xy
    else
        du_dy = 0.0
    end

    # calculating du_dtheta
    if obs_3p == false && obs_3n == false
        du_dtheta = (v_3p - v_3n)/(2*env.h_theta)
    elseif obs_3p == true && obs_3n == false
        du_dtheta = (v_x - v_3n)/env.h_theta
    elseif obs_3p == false && obs_3n == true
        du_dtheta = (v_3p - v_x)/env.h_theta
    else
        du_dtheta = 0.0
    end

    dV = [du_dx, du_dy, du_dtheta]

    # println("dV: ", dV)

    return dV
end

# may produce asymmetry due to interpolation index selection
function value_interp(x::Vector{Float64}, U::Array{Float64, 3}, O::Array{Bool, 3}, env::Environment)
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

    # adjust theta within bounds
    if x_itp[3] > pi
        x_itp[3] -= 2*pi
    elseif x_itp[3] < -pi
        x_itp[3] += 2*pi
    end

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

    for i in (i_0, i_1)         # SPEED: allocation (thought tuple fixed this?)
        for j in (j_0, j_1)         # SPEED: allocation
            for k in (k_0, k_1)     # SPEED: allocation (shows up at different times, may depend on k_0 and k_1)
                # println("[i,j,k]: ", [i,j,k], ", [x,y,th]: ", [env.x_grid[i], env.y_grid[j], env.theta_grid[k]], ", U[i,j,k] = ", U[i,j,k])

                # if any nodes used for interp are in obstacle, whole interpolation for "x_1p" is thrown out
                if O[i,j,k] == true
                    obs_itp = true
                    u_itp = 0.0

                    return u_itp, obs_itp
                end
            end
        end
    end

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

# SPEED: gets called 18 times for every gradient calculation
function find_idx(val, array)
    if val in array
        idx = indexin(val, array)[1]    # SPEED: allocation (think it creates an array sometimes -> bad)
                                        #   - 101.8 μs when array, 8.7 μs when no array (in practice will never exactly be on grid node)
    else
        idx = searchsortedfirst(array, val) - 1
    end

    return idx
end