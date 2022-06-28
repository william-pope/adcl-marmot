

using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")

# planner hierarchy
#   - HJB_planner()
#       - HJB_action()
#           - value_gradient()
#               - value_interp()

# NOTE: not actually used in operation
function HJB_planner(x_0, U, dt, plan_steps, A, O, env::Environment, veh::Vehicle)
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
        u_k = HJB_action(x_k, U, A, O, env, veh)  # SPEED: focus is here
        push!(u_path, u_k)

        # simulate forward one time step
        x_k1 = runge_kutta_4(x_k, u_k, dt, car_EoM, veh)    # needs K_sub
        push!(x_path, x_k1)

        x_k = deepcopy(x_k1)
    end

    # println("steps in HJB path: ", step)

    return x_path, u_path, step
end


# dt = 0.1 vs dt = 1.0, what issues will appear? hitting obstacles? 
# reactive controller will modify velocity on path, how will this affect path execution and safety?

# ISSUE: front/back chattering is still an issue to be fixed
#   - action computed from dot product of value gradient and f(x,u)
#   - gradient of value function computed using center difference method and interpolation
#   - is issue caused by being near an obstacle? CDM may be sampling points on obstacle
#       - huge value along axis that lands in obstacle
#   - theta gradient seems like biggest difference between fordward/backward steps
#       - may just be this spot
#   - signs on vdot flip between fordward/backward actions at different steps

#   - change is caused by dV[1]*xdot[1]
#       - x interp point lands in obstacle
#       - probably doesn't happen every time an obstacle is touched because numbers balance out
#       - dV[1] is pretty similar between steps, but small variations in cos(theta) due to theta cause big swing in dV[1]*xdot[1]
#           - same variations in sin(theta), but multiplied by normal scale dV (problem is still obstacle values)

#   - obstacles should not be included in gradient calculations, since their value is meaningless
#   - if adjacent point is in obstacle, either:
#       - leave axis out of dot product (probably not best)
#       - take single point gradient from noon-obstacle direction

#   - does issue only occur at corners or anytime near an obstacle?
#       - seems like only at corner

#   - fix will be in interp_value()
#   - zeroing out affected partial seems to work, but single FDM calc would be more accurate

# ISSUE: new problem in open field, check for bug somewhere
#   - starts when state clears the corner, no longer DQing any interp nodes on +x side
#       - although it looks like a node should still be in the obstacle (very close though)
#       - x partial is positive when should be negative -> hitting obstacle (why not detected?)
#       - is it rounding up the y value? but obviously hitting obstacle
#       - (!): issue where orientation allows side of car to overlap corner of obstacle (yes)

# working state:
# x: [-21.480752122416565, 20.97053490637511, 0.5470138792735877]
# x_itp: [-20.480752122416565, 20.97053490637511, 0.5470138792735877]
# [i,j,k]: [31, 71, 23], [x,y,th]: [-20.0, 20.0, 0.6981317007977319]
# above node in obstacle

# issue state:
# x: [-21.467849234711082, 20.978183753374577, 0.5232049851451851]
# x_itp: [-20.467849234711082, 20.978183753374577, 0.5232049851451851]
# [i,j,k]: [31, 71, 21], [x,y,th]: [-20.0, 20.0, 0.34906585039886595]   val ~= 300
# [i,j,k]: [31, 71, 22], [x,y,th]: [-20.0, 20.0, 0.5235987755982987]    val ~= 300

# fixed state:
# x: [-21.467849234711082, 20.978183753374577, 0.5232049851451851]
# x_itp: [-20.467849234711082, 20.978183753374577, 0.5232049851451851]
# [i,j,k]: [31, 71, 21], [x,y,th]: [-20.0, 20.0, 0.34906585039886595], U[i,j,k] = 1000.0
# above node in obstacle


# check that gradient working normally once around corner
#   - looks good to me


# NOTE: called directly by POMDP rollout simulator
#   - need to optimize this function and its children

# calculates optimal action as a function of state
function HJB_action(x, U, A, O, env::Environment, veh::Vehicle)
    # TO-DO: make sure this isn't causing any issues, conflicts with other boundary mechanisms
    # adjust position when near edge for boundary issues
    if abs(x[1] - env.x_grid[1]) < 2*env.h_xy
        x[1] = env.x_grid[3]
    elseif abs(x[1] - env.x_grid[end]) < 2*env.h_xy
        x[1] = env.x_grid[end-2]
    end

    if abs(x[2] - env.y_grid[1]) < 2*env.h_xy
        x[2] = env.y_grid[3]
    elseif abs(x[2] - env.y_grid[end]) < 2*env.h_xy
        x[2] = env.y_grid[end-2]
    end

    # println("\n--- ---\nx: ", x)
    
    # calculate gradient at current state
    dV = value_gradient(x, U, O, env)       # SPEED: big focus

    # println("dV: ", dV)

    # calculate optimal action
    vdot_min = Inf
    a_opt = A[1]
    for i in 1:size(A,1)
        xs = SVector{3, Float64}(x)     # TO-DO: clean up variable types
        xdot = car_EoM(xs, A[i], veh)

        vdot = dot(dV, xdot)

        # temporary
        # if i == 1 || i == 4
            # println("a: ", A[i], ", vdot: ", vdot)
            # println("dV[1] * xdot[1] = ", dV[1]*xdot[1])
            # println("dV[2] * xdot[2] = ", dV[2]*xdot[2])
            # println("dV[3] * xdot[3] = ", dV[3]*xdot[3])
            # println("")
        # end

        if vdot < vdot_min
            vdot_min = vdot
            a_opt = A[i]
        end
    end

    # println("a_opt: ", a_opt, ", vdot_min: ", vdot_min)
    return a_opt
end

# TO-DO: make CDM its own function
    #   - if any interps return obstacle, need to throw out sample point
    #   - this is only use of value_interp, could combine with CDM if it helps

# NOTE: may want to address bigger question of time step/action set (HJB vs tree search) in determining true value
#   - may want to choose best action by propagating each forward for Dt
#   - okay if resulting paths are only optimal as Dt->0

# ISSUE: dealing with obstacles during gradient approximation
#   - unsure what to do when one of interpolation nodes is in an obstacle
#   - depending on grid size, interp node could be somewhat far away from actual point of interest
#   - need to somehow calculate a gradient without using affected nodes
#   - gradient process is only a function of the state, not related to actions being tested
#   - obstacles are not a part of the flow field, treat like blank regions
#   - partials can be set to 0 if needed, just implies a flat gradient along that axis
#       - probably more accurate to switch to single point difference method, but will be a little harder

function value_gradient(x::Vector{Float64}, U::Array{Float64, 3}, O::Array{Bool, 3}, env::Environment)
    # interpolate value at surrounding points
    v_1p, obs_1p = value_interp(x + [env.h_xy, 0.0, 0.0], U, O, env)
    v_1n, obs_1n = value_interp(x - [env.h_xy, 0.0, 0.0], U, O, env)    # SPEED: braodcast.jl?
    v_2p, obs_2p = value_interp(x + [0.0, env.h_xy, 0.0], U, O, env)
    v_2n, obs_2n = value_interp(x - [0.0, env.h_xy, 0.0], U, O, env)    # SPEED: allocation?
    v_3p, obs_3p = value_interp(x + [0.0, 0.0, env.h_theta], U, O, env)
    v_3n, obs_3n = value_interp(x - [0.0, 0.0, env.h_theta], U, O, env)

    # NOTE: need to test this for oblique obstacle edges
    # approximate partial derivatives using centered difference scheme
    any((obs_1p, obs_1n)) ? du_dx = 0.0 : du_dx = (v_1p - v_1n)/(2*env.h_xy)
    any((obs_2p, obs_2n)) ? du_dy = 0.0 : du_dy = (v_2p - v_2n)/(2*env.h_xy)
    any((obs_3p, obs_3n)) ? du_dtheta = 0.0 : du_dtheta = (v_3p - v_3n)/(2*env.h_theta)

    dV = [du_dx, du_dy, du_dtheta]

    return dV
end

# SPEED: value_interp() seems most consequential

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