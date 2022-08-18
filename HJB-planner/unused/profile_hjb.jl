# Hamilton-Jacobi-Bellman demonstration

using Plots
using StaticArrays
using BSON: @save, @load
using BenchmarkTools
using ProfileView
include("HJB_functions.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   
unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)   
veh = marmot

# define workspace
W = [[-3.010 -7.023];
    [3.010 -7.023];
    [3.010 7.023];
    [-3.010 7.023]]

# W = [[-3.0 -3.0];
#     [3.0 -3.0];
#     [3.0 3.0];
#     [-3.0 3.0]]

# W = [[-4.0 -4.0];
#     [4.0 -4.0];
#     [4.0 4.0];
#     [-4.0 4.0]]

# define target set
T_xy = [[-0.5 4.5];
        [0.5 4.5];
        [0.5 5.5];
        [-0.5 5.5]]

# T_xy = [[-1.0 -1.0];
#         [1.0 -1.0];
#         [1.0 1.0];
#         [-1.0 1.0]]

# T_xy = [[-0.75 -0.75];
#         [0.75 -0.75];
#         [0.75 0.75];
#         [-0.75 0.75]]

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

OA = [[-1.0 -1.5];
    [1.0 -1.5];
    [1.0 1.5];
    [-1.0 1.5]]

O_vec = [O1, O2, O3]
# O_vec = [OA]

# initialize state grid
h_xy = 0.25
h_theta = deg2rad(10)

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

# ISSUE: does value estimate match actual path length (time)
#   - is this ok/expected? look into further

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/adcl/..."
algs_path = algs_path_mac

# calculate HJB
run_HJB = true

if run_HJB == true
    du_tol = 0.01
    max_steps = 5000
    anim_bool = false
    
    @btime solve_HJB_PDE(du_tol, max_steps, env, veh, anim_bool)
    ProfileView.@profview solve_HJB_PDE(du_tol, max_steps, env, veh, anim_bool)

    # @save algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
    # @save algs_path*"HJB-planner/bson/env.bson" env
    # @save algs_path*"HJB-planner/bson/veh.bson" veh
else
    @load algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
end

N_grid = size(env.x_grid,1)*size(env.y_grid,1)*size(env.theta_grid,1)
println("total grid nodes = ", N_grid)


# 4) PLOTS --- --- ---

# # plot U as heat map
# anim = @animate for k_plot in 1:size(env.theta_grid,1)
#     p_k = heatmap(env.x_grid, env.y_grid, 
#                 transpose(U_HJB[:,:,k_plot]), clim=(0,15),
#                 # xlim=(-3.5,5.5),
#                 aspect_ratio=:equal, 
#                 # size=(775,1050),
#                 size=(800,1050),
#                 xlabel="x-axis [m]", ylabel="y-axis [m]", 
#                 # title="HJB Value Function",
#                 titlefontsize = 20,
#                 colorbar_title = "time-to-target [s]",
#                 legend=false, colorbar=false,
#                 # legend=:topright,
#                 legend_font_pointsize = 11,
#                 top_margin = -30*Plots.mm,
#                 left_margin = -8*Plots.mm,
#                 bottom_margin = 8*Plots.mm)

#     # p_k = plot(xlim=(-3.5,5.5),
#     #     aspect_ratio=:equal, size=(750,1050),
#     #     xlabel="x-axis [m]", ylabel="y-axis [m]", 
#     #     titlefontsize = 20,
#     #     legend=:topright,
#     #     legend_font_pointsize = 11,
#     #     top_margin = -30*Plots.mm,
#     #     left_margin = 8*Plots.mm,
#     #     bottom_margin = 0*Plots.mm)

#     plot_polygon(p_k, env.W, 2, :black, "Workspace")
#     plot_polygon(p_k, env.T_xy, 2, :green, "Target Set")
#     plot_polygon(p_k, env.O_vec[1], 2, :red, "Obstacle")
#     plot_polygon(p_k, env.O_vec[2], 2, :red, "")
#     plot_polygon(p_k, env.O_vec[3], 2, :red, "")

#     # vehicle figure
#     x_pos = 4.5
#     y_pos = -1

#     x_max = x_pos + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
#     y_min = y_pos - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

#     y = [x_pos, y_pos, env.theta_grid[k_plot]]
    
#     V_c = pose_to_corners(y, veh)
#     V = [[V_c[1][1] V_c[1][2]];
#         [V_c[2][1] V_c[2][2]];
#         [V_c[3][1] V_c[3][2]];
#         [V_c[4][1] V_c[4][2]]]
        
#     plot!(p_k, [x_max], [y_pos], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
#     plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
#     plot_polygon(p_k, V, 2, :blue, "Vehicle Orientation")

#     theta_deg = round(rad2deg(y[3]), digits=1)
#     annotate!(x_pos, y_pos+2, text("theta [deg]:\n$theta_deg", 14))

#     display(p_k)
# end

# # gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)


# # plot optimal paths
# dt = 0.01

# p_path = plot(aspect_ratio=:equal, size=(750,1050), 
#             xlabel="x-axis [m]", ylabel="y-axis [m]",
#             xlim=(-3.5,5.5),
#             titlefontsize = 20,
#             legend_font_pointsize = 11,
#             legend=false,
#             top_margin = -4*Plots.mm,
#             left_margin = 8*Plots.mm)

# plot_polygon(p_path, env.W, 2, :black, "Workspace")
# plot_polygon(p_path, env.T_xy, 2, :green, "Target Set")
# plot_polygon(p_path, env.O_vec[1], 2, :red, "Obstacle")
# plot_polygon(p_path, env.O_vec[2], 2, :red, "")
# plot_polygon(p_path, env.O_vec[3], 2, :red, "")

# Y_0 = [[0.0, -6.0, 1/2*pi]]
#         # [-1.25, -1.0, -7/12*pi],
#         # [1.5, -4.5, 1/12*pi]]

# # iterate through given initial poses
# for y_0 in Y_0
#     @time y_path, u_path, step = HJB_planner(y_0, U_HJB, dt, env, veh)

#     path_time = step*dt
#     println("path execution time: ", path_time, " sec")

#     plot!(p_path, getindex.(y_path,1), getindex.(y_path,2),
#         linewidth = 2,
#         label="Optimal Path")

#     plot!(p_path, [y_path[1][1]], [y_path[1][2]], 
#         markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
#         label="")

#     V_c = pose_to_corners(y_0, veh)
#     V = [[V_c[1][1] V_c[1][2]];
#         [V_c[2][1] V_c[2][2]];
#         [V_c[3][1] V_c[3][2]];
#         [V_c[4][1] V_c[4][2]]]

#     plot_polygon(p_path, V, 2, :blue, "")

#     plot!(p_path, [y_path[end][1]], [y_path[end][2]], 
#         markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
#         label="")

#     V_c = pose_to_corners(y_path[end], veh)
#     V = [[V_c[1][1] V_c[1][2]];
#         [V_c[2][1] V_c[2][2]];
#         [V_c[3][1] V_c[3][2]];
#         [V_c[4][1] V_c[4][2]]]

#     plot_polygon(p_path, V, 2, :blue, "")
# end

# display(p_path)