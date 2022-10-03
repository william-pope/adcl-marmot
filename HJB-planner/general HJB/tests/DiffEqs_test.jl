# DiffEqs_test.jl

using DifferentialEquations
using StaticArrays
using BenchmarkTools

# actually slow as hell somehow...
#   - hates MVectors too lol

function bicycle_3d_EoM!(x_dot, x, u, t)
    x_dot[1] = u[1]*cos(x[3])
    x_dot[2] = u[1]*sin(x[3])
    x_dot[3] = u[1]*(1/0.5)*tan(u[2])
end

# fun = ODEFunction(bicycle_3d_EoM)

x_k = MVector(1.0, 2.0, pi/3)
u_k = SA[1.0, 0.475]

dt = 0.5
tspan = (0.0, dt)

prob = ODEProblem(bicycle_3d_EoM!, x_k, tspan, u_k)

@btime sol = solve(prob, Tsit5())

@btime x_k1 = last(solve(prob).u)