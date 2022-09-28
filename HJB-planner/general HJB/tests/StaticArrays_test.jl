# StaticArrays_test.jl

using StaticArrays
using BenchmarkTools

function EoM_old(x::Vector{Float64}, u::Vector{Float64}, wb)
    xdot = zeros(Float64, 3)
    
    xdot[1] = u[1]*cos(x[3])    # x
    xdot[2] = u[1]*sin(x[3])    # y
    xdot[3] = u[1]*(1/wb)*tan(u[2])  # theta

    return xdot
end

function EoM_static(x::SVector{3, Float64}, u::SVector{2, Float64}, wb::Float64)
    xdot1 = u[1]*cos(x[3])
    xdot2 = u[1]*sin(x[3])
    xdot3 = u[1]*(1/wb)*tan(u[2])

    n = 3
    xdot = SVector{n, Float64}(xdot1, xdot2, xdot3)

    return xdot
end


# xs = SVector{3, Float64}(1.0, 2.5, pi/3)

# @btime begin
#     x = zeros(Float64, 3)
#     x[1] = 1.0
#     x[2] = 2.5
#     x[3] = pi/3
# end

# @btime begin
#     xs = SVector{3, Float64}(1.0, 2.5, pi/3)
# end

x = [1.0, 2.5, pi/3]
u = [1.0, 0.475]

xs = SVector{3, Float64}(1.0, 2.5, pi/3)
us = SVector{2, Float64}(1.0, 0.475)

wb = 0.5

@btime xdot = EoM_old(x, u, wb)

@btime xdots = EoM_static(xs, us, wb)