# HJB_utils.jl

using LazySets

function common_prop_HJB(x_k, a_k, Dt, substeps)
    x_k1_subpath = Array{Array{Float64, 1}, 1}(undef, substeps)

    Dt_sub = Dt / substeps

    x_kk = x_k

    # step through substeps from x_k
    for kk in 1:substeps
        # Dv applied on first substep only
        kk == 1 ? a_kk = a_k : a_kk = [a_k[1], 0.0]
            
        # propagate for Dt_sub
        x_kk1 = common_4d_DT_EoM_HJB(x_kk, a_kk, Dt_sub)

        # store new point in subpath
        x_k1_subpath[kk] = x_kk1

        # pass state to next loop
        x_kk = x_kk1
    end

    x_k1 = x_k1_subpath[end]

    return x_k1, x_k1_subpath
end

function common_4d_DT_EoM_HJB(x_k, a_k, Dt)
    # system parameters
    l = 0.324

    # break out current state
    xp_k = x_k[1]
    yp_k = x_k[2]
    theta_k = x_k[3]
    v_k = x_k[4]

    # break out action
    phi_k = a_k[1]
    Dv_k = a_k[2]

    # calculate change in state over discrete time interval
    xp_dot_k = (v_k + Dv_k) * cos(theta_k)
    yp_dot_k = (v_k + Dv_k) * sin(theta_k)
    theta_dot_k = (v_k + Dv_k) * 1/l * tan(phi_k)

    # calculate next state
    xp_k1 = xp_k + (xp_dot_k * Dt)
    yp_k1 = yp_k + (yp_dot_k * Dt)
    theta_k1 = theta_k + (theta_dot_k * Dt)
    v_k1 = v_k + (Dv_k)

    # ISSUE: mod not working from neg to pos
    theta_k1 = theta_k1 % (2*pi)
    if theta_k1 > pi
        theta_k1 -= 2*pi
    elseif theta_k1 < -pi
        theta_k1 += 2*pi
    end

    # reassemble state vector
    x_k1 = SA[xp_k1, yp_k1, theta_k1, v_k1]

    return x_k1
end

function optimize_action(x, a_ind_array, actions, value_array, Dt, sg)
    val_x = Inf
    a_ind_opt = 1
    
    for a_ind in a_ind_array
        a = actions[a_ind]

        cost_x_a = get_cost(x, a, Dt)

        x_p, _ = common_prop_HJB(x, a, Dt, 4)
        val_xp = interp_state_value(x_p, value_array, sg)

        qval_x_a = cost_x_a + val_xp

        if qval_x_a < val_x
            val_x = qval_x_a
            a_ind_opt = a_ind
        end
    end

    return a_ind_opt, val_x
end

function interp_state_value(x, value_array, sg)
    # check if current state is within state space
    for d in eachindex(x)
        if x[d] < sg.state_grid.cutPoints[d][1] || x[d] > sg.state_grid.cutPoints[d][end]
            val_x = 1e5
            
            return val_x
        end
    end

    # interpolate value at given state
    val_x = interpolate(sg.state_grid, value_array, x)

    return val_x
end

# (?): wonder if this will have the same issue as setting hard constraints on obstacle nodes
#   - maybe just return worst value from subpath if it passes obstacle threshold?
#   - yep, looks like during solving everything gets set to 1e5
#   - try just returning 

# (!): function name is somewhat misleading, since it really returns a single value representative of the action
function interp_subpath_value(x, x_subpath, value_array, sg)
    # iterate through each subpoint up to x
    for x_sub in x_subpath[1:end-1]
        val_x_sub = interp_state_value(x_sub, value_array, sg)

        # if subpoint is in collision, return bad value for end state
        if val_x_sub > 50.0
            val_x = val_x_sub
            return val_x
        end
    end

    # if path is collision-free, return actual value for end state
    val_x = interp_state_value(x, value_array, sg)

    return val_x
end

# used for GridInterpolations.jl indexing
function multi2single_ind(ind_m, sg)
    ind_s = 1
    for d in eachindex(ind_m)
        ind_s += (ind_m[d]-1)*prod(sg.state_grid.cut_counts[1:(d-1)])
    end

    return ind_s
end

# workspace checker
function in_workspace(x, env, veh)
    veh_body = state_to_body(x, veh)

    if issubset(veh_body, env.workspace)
        return true
    end
        
    return false
end

# obstacle set checker
function in_obstacle_set(x, env, veh)
    veh_body = state_to_body(x, veh)
    
    for obstacle in env.obstacle_list
        if isempty(intersection(veh_body, obstacle)) == false
            return true
        end
    end

    return false
end

# target set checker
function in_target_set(x, env, veh)
    state = Singleton(x[1:2])

    if issubset(state, env.goal)
        return true
    end
    
    # veh_body = state_to_body(x, veh)

    # if issubset(veh_body, env.goal)
    #     return true
    # end

    return false
end

# vehicle body transformation function
function state_to_body(x, veh)
    # rotate body about origin by theta
    rot_matrix = [cos(x[3]) -sin(x[3]); sin(x[3]) cos(x[3])]
    body = linear_map(rot_matrix, veh.origin_body)

    # translate body from origin by [x, y]
    trans_vec = x[1:2]
    LazySets.translate!(body, trans_vec)

    return body
end

# used to create circles as polygons in LazySets.jl
function VPolyCircle(cent_cir, r_cir)
    # number of points used to discretize edge of circle
    pts = 12

    # circle radius is used as midpoint radius for polygon faces (over-approximation)
    r_poly = r_cir/cos(pi/pts)

    theta_rng = range(0, 2*pi, length=pts+1)

    cir_vertices = [[cent_cir[1] + r_poly*cos(theta), cent_cir[2] + r_poly*sin(theta)] for theta in theta_rng]
    
    poly_cir = VPolygon(cir_vertices)
    
    return poly_cir
end