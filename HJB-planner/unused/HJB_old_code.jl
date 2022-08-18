# HJB_old_code.jl

# finite difference value update
function update_node_value_FDM(value_array, init_array, i::Int, j::Int, k::Int, actions, env::Environment, veh::Vehicle, EoM::Function)   
    x1 = env.x_grid[i]
    x2 = env.y_grid[j]
    x3 = env.theta_grid[k]

    x = SVector{3, Float64}(x1, x2, x3)

    # compute val_ijk update
    val_ijk_min = Inf
    init_ijk = false

    # code in this loop gets iterated the most
    for u in actions
        xdot = EoM(x, u, veh)
    
        # calculate upwind indices
        i_uw = i + Int(sign(xdot[1]))
        j_uw = j + Int(sign(xdot[2]))
        k_uw = k + Int(sign(xdot[3]))

        k_uw == size(env.theta_grid,1)+1 ? k_uw = 2 : k_uw = k_uw
        k_uw == 0 ? k_uw = size(env.theta_grid,1)-1 : k_uw = k_uw

        if any((init_array[i_uw,j,k], init_array[i,j_uw,k], init_array[i,j,k_uw])) == true
            # pull value from upwind points
            val_i_uw = value_array[i_uw, j, k]
            val_j_uw = value_array[i, j_uw, k]
            val_k_uw = value_array[i, j, k_uw]

            # calculate value for given action
            val_ijk_u = finite_diff_eqn(xdot, val_i_uw, val_j_uw, val_k_uw, env)

            if val_ijk_u < val_ijk_min
                val_ijk_min = val_ijk_u
            end

            init_ijk = true
        end
    end

    if init_ijk == true
        val_ijk = val_ijk_min
    else
        val_ijk = value_array[i,j,k]
    end

    return val_ijk, init_ijk
end

function finite_diff_eqn(xdot::SVector{3, Float64}, val_i_uw::Float64, val_j_uw::Float64, val_k_uw::Float64, env::Environment)
    s1 = sign(xdot[1])
    s2 = sign(xdot[2])
    s3 = sign(xdot[3])

    num = 1.0 + s1/env.h_xy*xdot[1]*val_i_uw + s2/env.h_xy*xdot[2]*val_j_uw + s3/env.h_theta*xdot[3]*val_k_uw
    den = s1/env.h_xy*xdot[1] + s2/env.h_xy*xdot[2] + s3/env.h_theta*xdot[3]

    val_ijk = num/den

    return val_ijk
end

# Takei functions
# computes value update for grid point ijk using HJB finite difference scheme
function update_value_UCLA(U, i, j, k, env::Environment, veh::Vehicle)
    theta_k = env.theta_grid[k]

    xi_k = Int(sign(cos(theta_k)))
    nu_k = Int(sign(sin(theta_k)))

    ip = i + xi_k
    in = i - xi_k
    jp = j + nu_k
    jn = j - nu_k
    k == size(env.theta_grid,1)-1 ? kp = 1 : kp = k + 1
    k == 1 ? kn = size(env.theta_grid,1)-1 : kn = k - 1

    u_ip = U[ip,j,k] 
    u_in = U[in,j,k]
    u_jp = U[i,jp,k]
    u_jn = U[i,jn,k]
    u_kp = U[i,j,kp]
    u_kn = U[i,j,kn]

    # compute u_ijk update
    G_p1 = G_p1_ijk(theta_k, u_ip, u_jp, u_kp, u_kn, env, veh)
    F_p1 = F_p1_ijk(theta_k, u_ip, u_jp, env, veh)
    G_n1 = G_n1_ijk(theta_k, u_in, u_jn, u_kp, u_kn, env, veh)
    F_n1 = F_n1_ijk(theta_k, u_in, u_jn, env, veh)

    up_ijk = min(G_p1, F_p1, G_n1, F_n1, U[i,j,k])

    return up_ijk
end

# finite difference approximations
# (eq 30)
function G_p1_ijk(theta_k, u_ip, u_jp, u_kp, u_kn, env::Environment, veh::Vehicle)
    rho = veh.wb/tan(veh.c_phi)

    num = env.h_xy/veh.c_vf + abs(cos(theta_k))*u_ip + abs(sin(theta_k))*u_jp + env.h_xy/(rho*env.h_theta)*min(u_kp, u_kn)
    den = abs(cos(theta_k)) + abs(sin(theta_k)) + env.h_xy/(rho*env.h_theta)
    G_p1 = num/den

    return G_p1
end

# (eq 31)
function F_p1_ijk(theta_k, u_ip, u_jp, env::Environment, veh::Vehicle)
    num = env.h_xy/veh.c_vf + abs(cos(theta_k))*u_ip + abs(sin(theta_k))*u_jp
    den = abs(cos(theta_k)) + abs(sin(theta_k))
    F_p1 = num/den

    return F_p1
end

# (eq 32)
function G_n1_ijk(theta_k, u_in, u_jn, u_kp, u_kn, env::Environment, veh::Vehicle)
    rho = veh.wb/tan(veh.c_phi)

    num = env.h_xy/veh.c_vb + abs(cos(theta_k))*u_in + abs(sin(theta_k))*u_jn + env.h_xy/(rho*env.h_theta)*min(u_kp, u_kn)
    den = abs(cos(theta_k)) + abs(sin(theta_k)) + env.h_xy/(rho*env.h_theta)
    G_n1 = num/den

    return G_n1
end

# (eq 33)
function F_n1_ijk(theta_k, u_in, u_jn, env::Environment, veh::Vehicle)
    num = env.h_xy/veh.c_vb + abs(cos(theta_k))*u_in + abs(sin(theta_k))*u_jn
    den = abs(cos(theta_k)) + abs(sin(theta_k))
    F_n1 = num/den

    return F_n1
end