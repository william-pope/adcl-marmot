# Hamilton-Jacobi-Bellman demonstration

using Plots
using BSON: @save, @load
using BenchmarkTools
using ProfileView

include("HJB_phi_solver_functions.jl")
include("HJB_phi_planner_functions.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 4.553, 0.324, 0.5207, 0.2762, 0.0889)   
unit_car = Vehicle(1.0, 0.5, 0.5, 5.0, 0.5, 0.75, 0.375, 0.125)   
veh = marmot

Am = [[a_v, a_phi_dot] for a_v in [-veh.c_vb, veh.c_vf], a_phi_dot in [-veh.c_phi_dot, 0.0, veh.c_phi_dot]]
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
T_xy = [[2.35 10.2];
        [3.15 10.2];
        [3.15 11.0];
        [2.35 11.0]]

# T_xy = [[95.0 83.0];
#         [100.0 83.0];
#         [100.0 87.0];
#         [95.0 87.0]]

T_theta = [[-pi, pi]]

# define obstacles
# - circular obstacles defined as [x, y, r], converted to polygon overapproximation
OC1 = circle_to_polygon([1.5, 3.0, 0.6])
OC2 = circle_to_polygon([4.0, 5.0, 0.6])
OC3 = circle_to_polygon([2.8, 8.0, 0.6])

# OC1 = circle_to_polygon([65.0, 35.0, 20.0])
# OC2 = circle_to_polygon([35.0, 70.0, 10.0])

O_vec = [OC1, OC2, OC3]
# O_vec = []

# initialize state grid
h_xy = 0.1
h_theta = deg2rad(10)
h_phi = 2*veh.c_phi/(10)

env = Environment(h_xy, 
                h_theta,
                h_phi,
                minimum(W[:,1]) : h_xy : maximum(W[:,1]),
                minimum(W[:,2]) : h_xy : maximum(W[:,2]),
                -pi : h_theta : pi,
                -(veh.c_phi + h_phi) : h_phi : veh.c_phi + h_phi,
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
run_HJB = false
plot_U_HJB = false
plot_paths = true

if run_HJB == true
    du_tol = 0.01
    max_steps = 5000
    anim_bool = false
    # @btime solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)
    U_HJB, target_mat, obstacle_mat = solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)

    N_grid = size(env.x_grid,1) * size(env.y_grid,1) * size(env.theta_grid,1) * size(env.phi_grid,1)
    println("total grid nodes = ", N_grid)

    @save algs_path*"HJB-planner/bson/U_HJB_phi.bson" U_HJB
    @save algs_path*"HJB-planner/bson/target_mat_phi.bson" target_mat
    @save algs_path*"HJB-planner/bson/obstacle_mat_phi.bson" obstacle_mat
    @save algs_path*"HJB-planner/bson/env_phi.bson" env
    @save algs_path*"HJB-planner/bson/veh_phi.bson" veh
else
    @load algs_path*"HJB-planner/bson/U_HJB_phi.bson" U_HJB
    @load algs_path*"HJB-planner/bson/target_mat_phi.bson" target_mat
    @load algs_path*"HJB-planner/bson/obstacle_mat_phi.bson" obstacle_mat
end

# # @btime HJB_action(x_0, U_HJB, A, obstacle_mat, env, veh)

# # ProfileView.@profview for _ in 1:1000
# #     HJB_action(x_0, U_HJB, A, obstacle_mat, env, veh)
# # end

# generate optimal path
X_0 = [[2.75, 0.5, 0.5*pi, 0.0]]
# ,
#     [2.75, 0.5, 0.5*pi, veh.c_phi],
#     [2.75, 0.5, 0.5*pi, -veh.c_phi]]

dt = 0.01
plan_steps = 1000

# # path_time = step*dt
# # println("path execution time: ", path_time, " sec")


# 4) PLOTS --- --- ---

if plot_U_HJB == true
    # plot U as heat map
    # anim = @animate 

    i_plot = 29
    j_plot = 10

    for k_plot in [28] #1:size(env.theta_grid,1) #[28]
        println("at state: [$(env.x_grid[i_plot]), $(env.y_grid[j_plot]), $(env.theta_grid[k_plot]), :]")

        for p_plot in 2:size(env.phi_grid,1)-1 #[7]
            p_k = heatmap(env.x_grid, env.y_grid, 
                        transpose(U_HJB[:, :, k_plot, p_plot]), clim=(0,10),
                        aspect_ratio=:equal, 
                        size=(800,800),
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

            println("phi: ", env.phi_grid[p_plot], ",\t U: ", U_HJB[i_plot, j_plot, k_plot, p_plot])
            plot!(p_k, [env.x_grid[i_plot]], [env.y_grid[j_plot]], 
                markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

            # p_k = plot(aspect_ratio=:equal, 
            #             size=(750,1000),
            #             xlim=(3.2, 5.2),
            #             ylim=(0, 2.5),
            #             # xlabel="x-axis [m]", ylabel="y-axis [m]", 
            #             # title="HJB Value Function",
            #             titlefontsize = 20,
            #             # legend=:topright,
            #             legend=false, 
            #             legend_font_pointsize = 11,
            #             top_margin = -30*Plots.mm,
            #             left_margin = 8*Plots.mm,
            #             bottom_margin = 4*Plots.mm)

            # for i in 1:size(env.x_grid,1)
            #     for j in 1:size(env.y_grid,1)
            #         x_i = env.x_grid[i]
            #         y_j = env.y_grid[j]

            #         plot!(p_k, [x_i], [y_j],
            #             markershape=:circle, markersize=4, markercolor=:black,
            #             label="")
            #     end
            # end

            plot_polygon(p_k, env.W, 3, :black, "Workspace")
            plot_polygon(p_k, env.T_xy, 3, :green, "Target Set")
            # plot_polygon(p_k, env.O_vec[1], 3, :red, "Obstacle")
            for O in env.O_vec
                plot_polygon(p_k, O, 3, :red, "")
            end

            # vehicle figure
            x_pos = 6.75
            y_pos = 4.5

            x_max = x_pos + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
            y_min = y_pos - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

            x = [x_pos, y_pos, env.theta_grid[k_plot], env.phi_grid[p_plot]]
            
            V_c = pose_to_edges(x, veh)
            V = [[V_c[1][1] V_c[1][2]];
                [V_c[2][1] V_c[2][2]];
                [V_c[3][1] V_c[3][2]];
                [V_c[4][1] V_c[4][2]]]
                
            plot!(p_k, [x_pos], [y_pos], 
                markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot_polygon(p_k, V, 2, :blue, "Vehicle")

            theta_deg = round(rad2deg(x[3]), digits=1)
            annotate!(x_pos, y_pos+2.5, text("theta [deg]:\n$theta_deg", 14))

            phi_deg = round(rad2deg(x[4]), digits=1)
            annotate!(x_pos, y_pos+1.5, text("phi [deg]:\n$phi_deg", 14))

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
    # plot_polygon(p_path, env.O_vec[1], 3, :red, "Obstacle")
    for O in env.O_vec
        plot_polygon(p_path, O, 3, :red, "")
    end

    for x_0 in X_0
        # generate optimal path
        x_path, u_path, step = HJB_planner(x_0, U_HJB, dt, plan_steps, A, obstacle_mat, car_phi_EoM, env, veh)
        # display(u_path)

        p_phi = plot(getindex.(x_path,4),
            xlabel="time step", ylabel="phi [rad]",
            label="steering angle (phi)",
            legend=:top)

        plot!(p_phi, [0, 678], [veh.c_phi, veh.c_phi], linecolor=:black, linestyle=:dash, label="max steering angle")
        plot!(p_phi, [0, 678], [-veh.c_phi, -veh.c_phi], linecolor=:black, linestyle=:dash, label="")
        display(p_phi)

        path_time = step*dt
        println("path execution time: ", path_time, " sec")

        plot!(p_path, getindex.(x_path,1), getindex.(x_path,2),
            linewidth = 2, #linecolor=:purple,
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