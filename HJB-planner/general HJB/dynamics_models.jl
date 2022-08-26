# dynamics_models.jl

# equations of motion for 3 DoF kinematic bicycle model
# u = [u_v, u_phi]
function bicycle_3d_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 3)
    
    xdot[1] = u[1]*cos(x[3])    # x
    xdot[2] = u[1]*sin(x[3])    # y
    xdot[3] = u[1]*(1/veh.wheelbase)*tan(u[2])  # theta

    return xdot
end

# equations of motion for 4 DoF velocity kinematic bicycle model
# u = [u_a, u_phi]
function bicycle_4d_v_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 4)
    
    xdot[1] = x[4]*cos(x[3])    # x
    xdot[2] = x[4]*sin(x[3])    # y
    xdot[3] = x[4]*(1/veh.wheelbase)*tan(u[2])  # theta
    xdot[4] = u[1]  # v

    return xdot
end

# equations of motion for 4 DoF steering kinematic bicycle model
# u = [u_v, u_xi]
function bicycle_4d_phi_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 4)
    
    xdot[1] = u[1]*cos(x[3])    # x
    xdot[2] = u[1]*sin(x[3])    # y
    xdot[3] = u[1]*(1/veh.wheelbase)*tan(x[4])  # theta
    xdot[4] = u[2]  # phi

    return xdot
end

# equations of motion for 5 DoF kinematic bicycle model
# u = [u_a, u_xi]
function bicycle_5d_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 5)
    
    xdot[1] = x[4]*cos(x[3])    # x
    xdot[2] = x[4]*sin(x[3])    # y
    xdot[3] = x[4]*(1/veh.wheelbase)*tan(x[5])  # theta
    xdot[4] = u[1]  # v
    xdot[5] = u[2]  # phi

    return xdot
end