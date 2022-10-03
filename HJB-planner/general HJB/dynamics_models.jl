# dynamics_models.jl

using StaticArrays

# equations of motion for 3 DoF kinematic bicycle model
# x = [x, y, theta]; u = [u_v, u_phi]
function bicycle_3d_EoM(x::SVector{3, Float64}, u::SVector{2, Float64}, veh::Vehicle)
    x_dot1 = u[1]*cos(x[3])
    x_dot2 = u[1]*sin(x[3])
    x_dot3 = u[1]*(1/veh.wheelbase)*tan(u[2])

    x_dot = SVector(x_dot1, x_dot2, x_dot3)

    return x_dot
end

# equations of motion for 4 DoF velocity kinematic bicycle model
# x = [x, y, theta, v]; u = [u_a, u_phi]
function bicycle_4d_v_EoM(x::SVector{4, Float64}, u::SVector{2, Float64}, veh::Vehicle)    
    x_dot1 = x[4]*cos(x[3])
    x_dot2 = x[4]*sin(x[3])
    x_dot3 = x[4]*(1/veh.wheelbase)*tan(u[2])
    x_dot4 = u[1]

    x_dot = SVector(x_dot1, x_dot2, x_dot3, x_dot4)

    return x_dot
end

# equations of motion for 4 DoF steering kinematic bicycle model
# x = [x, y, theta, phi]; u = [u_v, u_xi]
function bicycle_4d_phi_EoM(x::SVector{4, Float64}, u::SVector{4, Float64}, veh::Vehicle)    
    x_dot1 = u[1]*cos(x[3])
    x_dot2 = u[1]*sin(x[3])
    x_dot3 = u[1]*(1/veh.wheelbase)*tan(x[4])
    x_dot4 = u[2]

    x_dot = SVector(x_dot1, x_dot2, x_dot3, x_dot4)

    return x_dot
end

# equations of motion for 5 DoF kinematic bicycle model
# u = [u_a, u_xi]
function bicycle_5d_EoM(x::SVector{5, Float64}, u::SVector{5, Float64}, veh::Vehicle)
    x_dot[1] = x[4]*cos(x[3])
    x_dot[2] = x[4]*sin(x[3])
    x_dot[3] = x[4]*(1/veh.wheelbase)*tan(x[5])
    x_dot[4] = u[1]
    x_dot[5] = u[2]

    x_dot = SVector(x_dot1, x_dot2, x_dot3, x_dot4)

    return x_dot
end

# # in-place version:
# function bicycle_3d_EoM!(x_dot::Vector{Float64}, x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
#     x_dot[1] = u[1]*cos(x[3])
#     x_dot[2] = u[1]*sin(x[3])
#     x_dot[3] = u[1]*(1/veh.wheelbase)*tan(u[2])
# end

# # equations of motion for 3 DoF kinematic bicycle model
# # u = [u_v, u_phi]
# function bicycle_3d_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
#     x_dot = zeros(Float64, 3)
    
#     x_dot[1] = u[1]*cos(x[3])    # x
#     x_dot[2] = u[1]*sin(x[3])    # y
#     x_dot[3] = u[1]*(1/veh.wheelbase)*tan(u[2])  # theta

#     return x_dot
# end