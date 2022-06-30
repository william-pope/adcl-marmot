# Hamilton-Jacobi-Bellman demonstration

using Plots
using BSON: @save, @load
using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_planner_functions.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   
unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)   
veh = marmot

Am = [[a_v,a_phi] for a_v in [-veh.c_vb, veh.c_vf], a_phi in [-veh.c_phi, 0.0, veh.c_phi]]
A = reshape(Am, (length(Am),1))
sort!(A, dims=1)

# define workspace
# W = [[-3.010 -7.023];
#     [3.010 -7.023];
#     [3.010 7.023];
#     [-3.010 7.023]]

# W = [[-3.0 -3.0];
#     [3.0 -3.0];
#     [3.0 3.0];
#     [-3.0 3.0]]

# W = [[-4.0 -4.0];
#     [4.0 -4.0];
#     [4.0 4.0];
#     [-4.0 4.0]]

W = [[0.0 0.0];
    [100.0 0.0];
    [100.0 100.0];
    [0.0 100.0]]

# define target set
# T_xy = [[-0.5 4.5];
#         [0.5 4.5];
#         [0.5 5.5];
#         [-0.5 5.5]]

# T_xy = [[-1.0 -1.0];
#         [1.0 -1.0];
#         [1.0 1.0];
#         [-1.0 1.0]]

# T_xy = [[-0.75 -0.75];
#         [0.75 -0.75];
#         [0.75 0.75];
#         [-0.75 0.75]]

T_xy = [[95.0 83.0];
        [100.0 83.0];
        [100.0 87.0];
        [95.0 87.0]]

T_theta = [[-pi, pi]]

# define obstacles
O1 = [[-2.0 2.0];
    [1.0 2.0];
    [1.0 2.8];
    [-2.0 2.8]]

O2 = [[1.2 -2.5];
    [2.0 -2.5];
    [2.0 0.5];
    [1.2 0.5]]

O3 = [[-1.5 -4.5];
    [-0.3 -4.5];
    [-0.3 -3.0];
    [-1.5 -3.0]]

OB = [[30.0 0.0];
    [100.0 0.0];
    [100.0 70.0];
    [30.0 70.0]]

# NEW: circular obstacles defined as [x, y, r], converted to polygon overapproximation
OC1 = circle_to_polygon([65.0, 35.0, 20.0])
OC2 = circle_to_polygon([35.0, 70.0, 10.0])

# O_vec = [O1, O2, O3]
O_vec = [OC1, OC2]

# initialize state grid
h_xy = 2.0
h_theta = deg2rad(15)

env = Environment(h_xy, 
                h_theta,
                minimum(W[:,1]) : h_xy : maximum(W[:,1]),
                minimum(W[:,2]) : h_xy : maximum(W[:,2]),
                -pi : h_theta : pi,
                W,
                T_xy,
                T_theta,
                O_vec)


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/adcl/..."
algs_path = algs_path_mac

# calculate HJB
run_HJB = true
plot_results = true

if run_HJB == true
    du_tol = 0.01
    max_steps = 5000
    anim_bool = false
    @btime U_HJB, T, O = solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)
    U_HJB, T, O = solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)

    N_grid = size(env.x_grid,1)*size(env.y_grid,1)*size(env.theta_grid,1)
    println("total grid nodes = ", N_grid)

    @save algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
    @save algs_path*"HJB-planner/bson/T.bson" T
    @save algs_path*"HJB-planner/bson/O.bson" O
    @save algs_path*"HJB-planner/bson/env.bson" env
    @save algs_path*"HJB-planner/bson/veh.bson" veh
else
    @load algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
    @load algs_path*"HJB-planner/bson/T.bson" T
    @load algs_path*"HJB-planner/bson/O.bson" O
end

# # generate optimal action for given state
# x_0 = [15.2, 5.1, 0.51*pi]

# @btime HJB_action(x_0, U_HJB, A, O, env, veh)

# ProfileView.@profview for _ in 1:1000
#     HJB_action(x_0, U_HJB, A, O, env, veh)
# end

# generate optimal path
x_0 = [15.2, 5.1, 0.51*pi]

dt = 0.01
plan_steps = 2e4

x_path, u_path, step = HJB_planner(x_0, U_HJB, dt, plan_steps, A, O, env, veh)
# ProfileView.@profview HJB_planner(x_0, U_HJB, dt, plan_steps, A, O, env, veh)

x_path, u_path = HJB_planner(x_0, U_HJB, dt, plan_steps, A, O, env, veh)

path_time = step*dt
println("path execution time: ", path_time, " sec")


# planner profiling
#   - initial btime for full planner: 159.059 ms (1,549,500 allocations)
#   - btime for HJB_action: 102.607 μs (273 allocations)
#   - need to reduce allocations


# 4) PLOTS --- --- ---

if plot_results == true
    # plot U as heat map
    anim = @animate for k_plot in 1:size(env.theta_grid,1)
        p_k = heatmap(env.x_grid, env.y_grid, 
                    transpose(U_HJB[:,:,k_plot]), clim=(0,100),
                    aspect_ratio=:equal, 
                    size=(1000,1100),
                    xlabel="x-axis [m]", ylabel="y-axis [m]", 
                    # title="HJB Value Function",
                    titlefontsize = 20,
                    colorbar_title = "time-to-target [s]",
                    legend=false, colorbar=false,
                    # legend=:topright,
                    legend_font_pointsize = 11,
                    top_margin = -30*Plots.mm,
                    left_margin = -8*Plots.mm,
                    bottom_margin = 8*Plots.mm)

        # p_k = plot(xlim=(-3.5,5.5),
        #     aspect_ratio=:equal, size=(750,1050),
        #     xlabel="x-axis [m]", ylabel="y-axis [m]", 
        #     titlefontsize = 20,
        #     legend=:topright,
        #     legend_font_pointsize = 11,
        #     top_margin = -30*Plots.mm,
        #     left_margin = 8*Plots.mm,
        #     bottom_margin = 0*Plots.mm)

        plot_polygon(p_k, env.W, 2, :black, "Workspace")
        plot_polygon(p_k, env.T_xy, 2, :green, "Target Set")
        for O in env.O_vec
            plot_polygon(p_k, O, 2, :red, "")
        end

        # vehicle figure
        x_pos = 110
        y_pos = 45

        x_max = x_pos + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
        y_min = y_pos - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

        x = [x_pos, y_pos, env.theta_grid[k_plot]]
        
        V_c = pose_to_edges(x, veh)
        V = [[V_c[1][1] V_c[1][2]];
            [V_c[2][1] V_c[2][2]];
            [V_c[3][1] V_c[3][2]];
            [V_c[4][1] V_c[4][2]]]
            
        plot!(p_k, [x_max], [y_pos], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        plot_polygon(p_k, V, 2, :blue, "Vehicle Orientation")

        theta_deg = round(rad2deg(x[3]), digits=1)
        annotate!(x_pos, y_pos+10, text("theta [deg]:\n$theta_deg", 14))

        display(p_k)
    end

    # gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)


    # plot path
    p_path = plot(aspect_ratio=:equal, size=(1400,1300), 
                xlabel="x-axis [m]", ylabel="y-axis [m]",
                titlefontsize = 20,
                legend_font_pointsize = 11,
                legend=false,
                top_margin = -4*Plots.mm,
                left_margin = 8*Plots.mm)

    plot_polygon(p_path, env.W, 2, :black, "Workspace")
    plot_polygon(p_path, env.T_xy, 2, :green, "Target Set")
    for O in env.O_vec
        plot_polygon(p_path, O, 2, :red, "")
    end

    plot!(p_path, getindex.(x_path,1), getindex.(x_path,2),
        linewidth = 2,
        label="Optimal Path")

    plot!(p_path, [x_path[1][1]], [x_path[1][2]], 
        markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
        label="")

    V_c = pose_to_edges(x_0, veh)
    V = [[V_c[1][1] V_c[1][2]];
        [V_c[2][1] V_c[2][2]];
        [V_c[3][1] V_c[3][2]];
        [V_c[4][1] V_c[4][2]]]

    plot_polygon(p_path, V, 2, :blue, "")

    plot!(p_path, [x_path[end][1]], [x_path[end][2]], 
        markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
        label="")

    V_c = pose_to_edges(x_path[end], veh)
    V = [[V_c[1][1] V_c[1][2]];
        [V_c[2][1] V_c[2][2]];
        [V_c[3][1] V_c[3][2]];
        [V_c[4][1] V_c[4][2]]]

    plot_polygon(p_path, V, 2, :blue, "")

    display(p_path)
end