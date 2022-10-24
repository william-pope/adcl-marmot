# utils.jl

using LazySets

struct Vehicle
    wheelbase::Float64
    body_dims::Array{Float64}
    origin_to_cent::Array{Float64}
    origin_body::VPolygon
end

# defines vehicle geometry
function define_vehicle(wheelbase, body_dims, origin_to_cent)
    x0_min = origin_to_cent[1] - 1/2*body_dims[1]
    x0_max = origin_to_cent[1] + 1/2*body_dims[1]
    y0_min = origin_to_cent[2] - 1/2*body_dims[2]
    y0_max = origin_to_cent[2] + 1/2*body_dims[2]
    origin_body = VPolygon([[x0_min, y0_min], [x0_max, y0_min], [x0_max, y0_max], [x0_min, y0_max]])

    veh = Vehicle(wheelbase, body_dims, origin_to_cent, origin_body)
    return veh
end

# equations of motion for 4 DoF velocity kinematic bicycle model
# x = [x, y, theta, v]; u = [u_a, u_phi]
function bicycle_4d_v_EoM(x, u, veh)    
    x_dot1 = x[4]*cos(x[3])
    x_dot2 = x[4]*sin(x[3])
    x_dot3 = x[4]*(1/veh.wheelbase)*tan(u[2])
    x_dot4 = u[1]

    x_dot = [x_dot1, x_dot2, x_dot3, x_dot4]

    return x_dot
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k, u_k, dt::Float64, EoM, veh)    
    w1 = EoM(x_k, u_k, veh)
    w2 = EoM(x_k + w1*dt/2, u_k, veh)
    w3 = EoM(x_k + w2*dt/2, u_k, veh)
    w4 = EoM(x_k + w3*dt, u_k, veh)

    x_k1 = x_k + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    return x_k1
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
    pts = 8

    # circle radius is used as midpoint radius for polygon faces (over-approximation)
    r_poly = r_cir/cos(pi/pts)

    theta_rng = range(0, 2*pi, length=pts+1)

    cir_vertices = [[cent_cir[1] + r_poly*cos(theta), cent_cir[2] + r_poly*sin(theta)] for theta in theta_rng]
    
    poly_cir = VPolygon(cir_vertices)
    
    return poly_cir
end