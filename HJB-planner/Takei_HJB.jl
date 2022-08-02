
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