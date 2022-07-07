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
W = [[0.0 0.0];
    [5.518 0.0];
    [5.518 11.036];
    [0.0 11.036]]

# W = [[0.0 0.0];
#     [100.0 0.0];
#     [100.0 100.0];
#     [0.0 100.0]]

# define target set
T_xy = [[2.25 10.0];
        [3.25 10.0];
        [3.25 11.0];
        [2.25 11.0]]

# T_xy = [[95.0 83.0];
#         [100.0 83.0];
#         [100.0 87.0];
#         [95.0 87.0]]

T_theta = [[-pi, pi]]

# define obstacles
# - circular obstacles defined as [x, y, r], converted to polygon overapproximation
OC1 = circle_to_polygon([1.5, 3.0, 0.5])
OC2 = circle_to_polygon([4.0, 6.0, 0.5])
OC3 = circle_to_polygon([2.2, 7.25, 0.5])

# OC1 = circle_to_polygon([65.0, 35.0, 20.0])
# OC2 = circle_to_polygon([35.0, 70.0, 10.0])

O_vec = [OC1, OC2, OC3]

# initialize state grid
h_xy = 0.5
h_theta = deg2rad(9)

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
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"

algs_path = algs_path_mac

# calculate HJB
run_HJB = true
plot_U_HJB = true
plot_paths = true

if run_HJB == true
    du_tol = 0.01
    max_steps = 5000
    anim_bool = false
    @btime solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)
    U_HJB, target_mat, obstacle_mat = solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)

    N_grid = size(env.x_grid,1)*size(env.y_grid,1)*size(env.theta_grid,1)
    println("total grid nodes = ", N_grid)

    # println(U_HJB[6, 12, 5])
    # println(U_HJB[28, 43, 13])
    # println(U_HJB[45, 92, 28])
   
    # original: 225 steps
    #   - 7.719296560087587
    #   - 5.677956746308499
    #   - 2.05270933919414

    # Gauss-Seidel: 154 steps
    #   - 7.710956256430729
    #   - 5.677604833255936
    #   - 2.052709011880812

    # actual Gauss-Seidel: 50 steps
    #   - 7.715912440064768
    #   - 5.681829778931739
    #   - 2.0529986173047634

    @save algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
    @save algs_path*"HJB-planner/bson/target_mat.bson" target_mat
    @save algs_path*"HJB-planner/bson/obstacle_mat.bson" obstacle_mat
    @save algs_path*"HJB-planner/bson/env.bson" env
    @save algs_path*"HJB-planner/bson/veh.bson" veh
else
    @load algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
    @load algs_path*"HJB-planner/bson/target_mat.bson" target_mat
    @load algs_path*"HJB-planner/bson/obstacle_mat.bson" obstacle_mat
end

# # generate optimal action for given state
# x_0 = [15.2, 5.1, 0.51*pi]

# @btime HJB_action(x_0, U_HJB, A, obstacle_mat, env, veh)

# ProfileView.@profview for _ in 1:1000
#     HJB_action(x_0, U_HJB, A, obstacle_mat, env, veh)
# end

# # generate optimal path
# x_0 = [5.0, 0.5, 0.7*pi]

X_0 = [[2.75, 0.5, 0.5*pi]]
# ,
#         [4.75, 3.0, deg2rad(105)],
#         [2.75, 6.2, deg2rad(-85)],
#         [0.7, 9.0, deg2rad(45)]]

dt = 0.01
plan_steps = 1000

# x_path, u_path, step = HJB_planner(x_0, U_HJB, dt, plan_steps, A, obstacle_mat, env, veh)
# # ProfileView.@profview HJB_planner(x_0, U_HJB, dt, plan_steps, A, obstacle_mat, env, veh)

# path_time = step*dt
# println("path execution time: ", path_time, " sec")


# planner profiling
#   - initial btime for full planner: 159.059 ms (1,549,500 allocations)
#   - btime for HJB_action: 102.607 μs (273 allocations)
#   - need to reduce allocations


# 4) PLOTS --- --- ---

if plot_U_HJB == true
    # plot U as heat map
    anim = @animate for k_plot in [31] #1:size(env.theta_grid,1)
        p_k = heatmap(env.x_grid, env.y_grid, 
                    transpose(U_HJB[:,:,k_plot]), clim=(0,10),
                    aspect_ratio=:equal, 
                    size=(600,1000),
                    # xlabel="x-axis [m]", ylabel="y-axis [m]", 
                    # title="HJB Value Function",
                    titlefontsize = 20,
                    # legend=:topright,
                    legend=false, 
                    colorbar_title = "time-to-target [s]",
                    colorbar=false,
                    legend_font_pointsize = 11,
                    top_margin = -30*Plots.mm,
                    left_margin = 8*Plots.mm,
                    bottom_margin = 4*Plots.mm)

        # p_k = plot(aspect_ratio=:equal, 
        #             size=(750,1000),
        #             # xlabel="x-axis [m]", ylabel="y-axis [m]", 
        #             # title="HJB Value Function",
        #             titlefontsize = 20,
        #             legend=:topright,
        #             # legend=false, 
        #             legend_font_pointsize = 11,
        #             top_margin = -30*Plots.mm,
        #             left_margin = 8*Plots.mm,
        #             bottom_margin = 4*Plots.mm)

        # for i in 1:size(env.x_grid,1)
        #     for j in 1:size(env.y_grid,1)
        #         x_i = env.x_grid[i]
        #         y_j = env.y_grid[j]

        #         plot!(p_k, [x_i], [y_j],
        #             markershape=:circle, markersize=2, markercolor=:black,
        #             label="")
        #     end
        # end

        plot_polygon(p_k, env.W, 3, :black, "Workspace")
        plot_polygon(p_k, env.T_xy, 3, :green, "Target Set")
        plot_polygon(p_k, env.O_vec[1], 3, :red, "Obstacle")
        for O in env.O_vec
            plot_polygon(p_k, O, 3, :red, "")
        end

        # # vehicle figure
        # x_pos = 6.75
        # y_pos = 5

        # x_max = x_pos + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
        # y_min = y_pos - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

        # x = [x_pos, y_pos, env.theta_grid[k_plot]]
        
        # V_c = pose_to_edges(x, veh)
        # V = [[V_c[1][1] V_c[1][2]];
        #     [V_c[2][1] V_c[2][2]];
        #     [V_c[3][1] V_c[3][2]];
        #     [V_c[4][1] V_c[4][2]]]
            
        # plot!(p_k, [x_max], [y_pos], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        # plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        # plot_polygon(p_k, V, 2, :blue, "Vehicle")

        # theta_deg = round(rad2deg(x[3]), digits=1)
        # annotate!(x_pos, y_pos+1.5, text("theta [deg]:\n$theta_deg", 14))

        # # adding path to heatmap
        # dt = 0.01
        # plan_steps = 1000

        # x_0 = [2.75, 0.5, 0.5*pi]
        # x_path, u_path, step = HJB_planner(x_0, U_HJB, dt, plan_steps, A, obstacle_mat, env, veh)

        # plot!(p_k, getindex.(x_path,1), getindex.(x_path,2),
        #     linewidth = 2, linecolor=:white,
        #     label="HJB Path")

        # plot!(p_k, [x_path[1][1]], [x_path[1][2]], 
        #     markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, 
        #     label="")

        # V_c = pose_to_edges(x_0, veh)
        # V = [[V_c[1][1] V_c[1][2]];
        #     [V_c[2][1] V_c[2][2]];
        #     [V_c[3][1] V_c[3][2]];
        #     [V_c[4][1] V_c[4][2]]]

        # plot_polygon(p_k, V, 2, :white, "Vehicle")

        # plot!(p_k, [x_path[end][1]], [x_path[end][2]], 
        #     markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, 
        #     label="")

        # V_c = pose_to_edges(x_path[end], veh)
        # V = [[V_c[1][1] V_c[1][2]];
        #     [V_c[2][1] V_c[2][2]];
        #     [V_c[3][1] V_c[3][2]];
        #     [V_c[4][1] V_c[4][2]]]

        # plot_polygon(p_k, V, 2, :white, "")

        display(p_k)
    end

    # gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)
end

if plot_paths == true
    # plot path
    p_path = plot(aspect_ratio=:equal, size=(600,1000),
                # xlabel="x-axis [m]", ylabel="y-axis [m]",
                titlefontsize = 20,
                legend_font_pointsize = 11,
                # legend=:topright,
                legend=false,
                top_margin = -30*Plots.mm,
                left_margin = 8*Plots.mm,
                bottom_margin = 4*Plots.mm)

    plot_polygon(p_path, env.W, 3, :black, "Workspace")
    plot_polygon(p_path, env.T_xy, 3, :green, "Target Set")
    plot_polygon(p_path, env.O_vec[1], 3, :red, "Obstacle")
    for O in env.O_vec
        plot_polygon(p_path, O, 3, :red, "")
    end

    for x_0 in X_0
        # generate optimal path
        x_path, u_path, step = HJB_planner(x_0, U_HJB, dt, plan_steps, A, obstacle_mat, env, veh)

        path_time = step*dt
        println("path execution time: ", path_time, " sec")

        plot!(p_path, getindex.(x_path,1), getindex.(x_path,2),
            linewidth = 2, linecolor=:purple,
            label="")

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
end