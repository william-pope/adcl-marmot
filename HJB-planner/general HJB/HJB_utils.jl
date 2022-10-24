# HJB_utils.jl

using LazySets

function common_prop_HJB(x_k, a_k, Dt, substeps)
    x_k1_subpath = Array{Array{Float64, 1}, 1}(undef, substeps)

    Dt_sub = Dt / substeps

    x_kk = x_k

    # step through substeps from x_k
    for kk in 1:substeps
        # Dv applied on the first substep only
        kk == 1 ? a_kk = a_k : a_kk = [0.0, a_k[2]]
            
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
    Dv_k = a_k[1]
    phi_k = a_k[2]

    # calculate change in state over discrete time interval
    xp_dot_k = (v_k + Dv_k) * cos(theta_k)
    yp_dot_k = (v_k + Dv_k) * sin(theta_k)
    theta_dot_k = (v_k + Dv_k) * 1/l * tan(phi_k)

    # calculate next state
    xp_k1 = xp_k + (xp_dot_k * Dt)
    yp_k1 = yp_k + (yp_dot_k * Dt)
    theta_k1 = theta_k + (theta_dot_k * Dt)
    v_k1 = v_k + (Dv_k)

    theta_k1 = theta_k1 % (2*pi)
    if theta_k1 > pi
        theta_k1 -= 2*pi
    end

    # reassemble state vector
    x_k1 = SA[xp_k1, yp_k1, theta_k1, v_k1]

    return x_k1
end

function interp_value(x, value_array, sg)
    # check if current state is within state space
    for d in eachindex(x)
        if x[d] < sg.state_grid.cutPoints[d][1] || x[d] > sg.state_grid.cutPoints[d][end]
            val_itp = 1e5
            return val_itp
        end
    end

    # interpolate value at given state
    val_itp = interpolate(sg.state_grid, value_array, x)

    return val_itp
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