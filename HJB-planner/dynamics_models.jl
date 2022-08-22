# dynamics_models.jl

using StaticArrays

struct Environment
    h_xy::Float64
    h_theta::Float64
    x_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    y_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    theta_grid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    W::Matrix{Float64}
    T_xy::Matrix{Float64}
    T_theta::Vector{Vector{Float64}}
    O_vec::Vector{Matrix{Float64}}
end

struct Vehicle
    u_vf_max::Float64   # max forward speed [m/s]
    u_vb_max::Float64   # max backward speed [m/s]
    u_phi_max::Float64  # max steering angle [rad]
    axle_l::Float64     # wheelbase [m]
    ext_l::Float64      # length [m]
    ext_w::Float64      # width [m]
    ext2axle::Float64    # rear bumber to rear axle [m]
end

# equations of motion for 3 DoF kinematic bicycle model
function dubins_car_EoM(x::SVector{3, Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot1 = u[1]*cos(x[3])
    xdot2 = u[1]*sin(x[3])
    xdot3 = u[1]*(1/veh.axle_l)*tan(u[2])
    xdot = SVector{3, Float64}(xdot1, xdot2, xdot3)
    return xdot
end
