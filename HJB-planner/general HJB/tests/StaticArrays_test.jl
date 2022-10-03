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

    xdot = SVector(xdot1, xdot2, xdot3)

    return xdot
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k, u_k, dt::Float64, EoM, veh, sg)    
    w1 = EoM(x_k, u_k, veh)
    w2 = EoM(x_k + w1*dt/2, u_k, veh)
    w3 = EoM(x_k + w2*dt/2, u_k, veh)
    w4 = EoM(x_k + w3*dt, u_k, veh)

    x_k1 = x_k + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    # x_k1 = mod_state_angle(x_k1, sg)

    return x_k1
end

# adjust angles within [-pi,pi] bounds
function mod_state_angle(x, sg)
    for d in eachindex(x)
        if sg.angle_wrap_array[d] == true
            x[d] = x[d] % (2*pi)
            x[d] > pi ? x[d] -= 2*pi : x[d] -= 0
        end
    end

    return x
end


function mod_angle(x::Float64, wrap::Bool)
    if wrap == true
        xm = x % (2*pi)
        xm > pi ? xm -= 2*pi : xm -= 0.0
    else
        xm = x
    end

    return xm
end

# xs = SVector(1.0, 2.5, 7/4*pi)

# wrap_array = [false, false, true]

# n = 3

# @btime (1,2,3)
# @btime collect(1:n)
# @btime tuple(1:$n...)

# xm = map(i -> mod_angle(xs[i], wrap_array[i]), SA[1,2,3])

xa = [[1.0, 2.5, 3/4*pi]]

@btime SA[1.0, 2.5, 3/4*pi]

@btime SA[$getindex($xa, $1)...]

@btime @SVector [getindex($xa, $1)...]

@btime SA[xa[1][:]]



# println("\nindexing StaticArray")
# @btime s = xs[1]

# println("\nindexing standard array")
# @btime a = xa[1]

# println("\nfunction with StaticArray")
# @btime map(i -> xs[i]^2, (1,2,3))

# println("\nfunction with standard array")
# @btime map(i -> xa[i]^2, (1,2,3))

# println("\nfunction without indexing")
# @btime map(i -> i^2, (1,2,3))


# @btime begin
#     x = zeros(Float64, 3)
#     x[1] = 1.0
#     x[2] = 2.5
#     x[3] = pi/3
# end

# @btime begin
#     xs = SVector{3, Float64}(1.0, 2.5, pi/3)
# end

# x = [1.0, 2.5, pi/3]
# u = [1.0, 0.475]

# xs = SVector{3, Float64}(1.0, 2.5, pi/3)
# us = SVector{2, Float64}(1.0, 0.475)

# wb = 0.5

# @btime xdot = EoM_old(x, u, wb)

# @btime xdots = EoM_static(xs, us, wb)