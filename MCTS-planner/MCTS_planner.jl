# final project

using Plots
using JLD

include("MCTS_functions.jl")
include("../HJB-planner/HJB_functions.jl")

# need for MDP:
#   - state space (S)
#   - action space (A)
#   - transition function (T(sp|s,a))
#   - reward function (R(s,a))
#   - discount (gamma)s


# 2) PARAMETERS --- --- ---

# vehicle parameters
unit_car = Vehicle(1.0, 0.75, 0.5, 0.5, 0.75, 0.25, 0.125)   

# define workspace
W_set = [[-10.0 -10.0];
        [10.0 -10.0];
        [10.0 10.0];
        [-10.0 10.0]]

# define target set
T_xy_set = [[-1.0 4.0];
            [1.0 4.0];
            [1.0 6.0];
            [-1.0 6.0]]

T_theta_set = [[-pi, pi]]

# define obstacles
O1_set = [[-2.0 -1.0];
        [3.5 -1.0];
        [3.5 1.0];
        [-2.0 1.0]]

O2_set = [[-6.5 -5.0];
        [-2.0 -5.0];
        [-2.0 -3.0];
        [-6.5 -3.0]]

O_set = [O1_set, O2_set]
# O_set = []

# initialize state grid
h_xy = 0.125
h_theta = deg2rad(2.5)

sg = StateGrid(h_xy, 
                h_theta,
                minimum(W_set[:,1]) : h_xy : maximum(W_set[:,1]),
                minimum(W_set[:,2]) : h_xy : maximum(W_set[:,2]),
                -pi : h_theta : pi)


# 3) MAIN --- --- ---
println("\n--- START ---")

# HJB
run_HJB = true

if run_HJB == true
    du_tol = 0.01
    max_reps = 100
    @time U = solve_HJB_PDE(du_tol, max_reps, sg, O_set, unit_car)
    V_HJB = -U

    n_nodes = size(sg.x_grid, 1) * size(sg.y_grid, 1) * size(sg.theta_grid, 1)
    println("total number of grid nodes: ", n_nodes)

    save("V_HJB.jld", "V_HJB", V_HJB)
else
    V_HJB = load("V_HJB.jld", "V_HJB")
    U = -V_HJB
end

@load "../HJB-planner/bson/U_HJB.bson" U_HJB
V_HJB = -U_HJB

@load "../HJB-planner/bson/env.bson" env
@load "../HJB-planner/bson/veh.bson" veh

# MCTS
S = [[-10.0, 10.0],
    [-10.0, 10.0],
    [-pi, pi]]

A_v = [-0.751, 1.0]
A_phi = [-0.5, 0.0, 0.5]
A = vec([[a_v, a_phi] for a_phi in A_phi, a_v in A_v])

gamma = 0.95

# ISSUE: noise works, but should add to action instead of x_dot
#   - needs different scaling for x,y (orientation) and theta (units)
#   - should be a "correct" way to add noise to bicycle model

# not doing anything?

s_0_mcts = [-4, -8, -0.55*pi]
dt = 0.1

println("planning path 1")
std_v = 0.02
std_phi = 0.01
sims = 324
w_max = 3

s_0_mcts = [-5, -7, -0.5*pi]

s_hist1_mcts, a_hist1_mcts = mcts_planner(s_0_mcts, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, car_EoM, sg, unit_car)

println("planning path 2")
std_v = 0.1
std_phi = 0.05
sims = 324
w_max = 3

s_0_mcts = [-3.75, -8.25, -0.5*pi]

s_hist2_mcts, a_hist2_mcts = mcts_planner(s_0_mcts, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, car_EoM, sg, unit_car)

println("planning path 3")
std_v = 0.5
std_phi = 0.25
sims = 324
w_max = 3

s_0_mcts = [-2.5, -9.5, -0.5*pi]

s_hist3_mcts, a_hist3_mcts = mcts_planner(s_0_mcts, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, car_EoM, sg, unit_car)

# s_0_hjb = [-4.0, -8, pi/2]
# y_path_HJB, u_path_HJB = HJB_planner(s_0_hjb, U, T_xy_set, T_theta_set, dt, sg, unit_car)

# display(s_hist_mcts)
# display(a_hist_mcts)


# 4) PLOTS --- --- ---

# plot U as heat map
if run_HJB == true
    for k_plot in LinRange(1, size(sg.theta_grid, 1), 15)
        k_plot = Int(round(k_plot, digits=0))
        theta_k = round(rad2deg(sg.theta_grid[k_plot]), digits=3)

        p_k = heatmap(sg.x_grid, sg.y_grid, transpose(V_HJB[:,:,k_plot]), clim=(-20,0),
                    aspect_ratio=:equal, size=(675,600),
                    xlabel="x-axis [m]", ylabel="y-axis [m]",
                    right_margin = 4Plots.mm,
                    top_margin = -8Plots.mm,
                    bottom_margin = -8Plots.mm)
                    # , title="HJB Value Function: u(x, y, theta=$theta_k)"

        plot_polygon(W_set, p_k, 3, :black, "Workspace")
        plot_polygon(T_xy_set, p_k, 3, :green, "Target Set")
        plot_polygon(O1_set, p_k, 3, :red, "Obstacle")
        plot_polygon(O2_set, p_k, 3, :red, "")

        display(p_k)
    end
end

# plot optimal path from y_0 to target set
p_path_mcts = plot(aspect_ratio=:equal, size=(600,600), legend=:topright)
            # xlabel="x-axis [m]", ylabel="y-axis [m]")

plot_polygon(W_set, p_path_mcts, 3, :black, "Workspace")
plot_polygon(T_xy_set, p_path_mcts, 3, :green, "Target Set")
plot_polygon(O1_set, p_path_mcts, 3, :red, "Obstacle")
plot_polygon(O2_set, p_path_mcts, 3, :red, "")

# plot!(p_path_mcts, getindex.(y_path_HJB,1), getindex.(y_path_HJB,2), 
#     color=:black, linewidth=2.5, linestyle=:dash, label="HJB Path")

plot!(p_path_mcts, getindex.(s_hist1_mcts,1), getindex.(s_hist1_mcts,2),
    color=:green, linewidth=2.5, label="MCTS, σ1")

plot!(p_path_mcts, getindex.(s_hist2_mcts,1), getindex.(s_hist2_mcts,2),
    color=:blue, linewidth=2.5, label="MCTS, σ2")

plot!(p_path_mcts, getindex.(s_hist3_mcts,1), getindex.(s_hist3_mcts,2),
    color=:purple, linewidth=2.5, label="MCTS, σ3")

display(p_path_mcts)