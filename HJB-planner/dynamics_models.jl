# dynamics_models.jl

using StaticArrays

# equations of motion for 3 DoF kinematic bicycle model
function dubins_car_EoM(x::SVector{3, Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot1 = u[1]*cos(x[3])
    xdot2 = u[1]*sin(x[3])
    xdot3 = u[1]*(1/veh.wb)*tan(u[2])

    xdot = SVector{3, Float64}(xdot1, xdot2, xdot3)

    return xdot
end