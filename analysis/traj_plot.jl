# trajectory plot for physical system

using Plots
using BSON: @save, @load
using CSV
using DataFrames
include("../HJB-planner/HJB_functions.jl")

hjb_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/HJB-planner/bson/"
@load hjb_path_mac*"env.bson" env
@load hjb_path_mac*"veh.bson" veh

hist_path_mac = "/Users/willpope/Desktop/Research/marmot-ros/controller_pkg/histories/stored/"
@load hist_path_mac*"s_hist.bson" s_hist

vrpn_csv = CSV.File(hist_path_mac*"s_hist_vrpn.csv")
s_hist_vrpn = []
for s_kv in vrpn_csv
    push!(s_hist_vrpn, s_kv)
end

@show size(s_hist,1)
@show size(s_hist_vrpn,1)

# plot optimal path from y_0 to target set
p_path = plot(aspect_ratio=:equal, size=(750,1050), 
            xlabel="x-axis [m]", ylabel="y-axis [m]",
            title="Vicon Output", 
            titlefontsize = 20,
            legend_font_pointsize = 11,
            top_margin = -5*Plots.mm,
            left_margin = 4*Plots.mm)

plot_polygon(p_path, env.W, 2, :black, "Workspace")
plot_polygon(p_path, env.T_xy, 2, :green, "Target Set")
# plot_polygon(p_path, env.O_vec[1], 2, :red, "Obstacle")
# plot_polygon(p_path, env.O_vec[2], 2, :red, "")
# plot_polygon(p_path, env.O_vec[3], 2, :red, "")]

plot!(p_path, getindex.(s_hist,1), getindex.(s_hist,2),
    linewidth=0, markershape=:circle, markersize=3, markerstrokewidth=0 ,
    label="Planning Points")

plot!(p_path, getindex.(s_hist_vrpn,1), getindex.(s_hist_vrpn,2),
    linewidth=1,
    label="Vicon Output")

# plot initial state
V_ci = pose_to_corners(s_hist[1], veh)
V_i = [[V_ci[1][1] V_ci[1][2]];
    [V_ci[2][1] V_ci[2][2]];
    [V_ci[3][1] V_ci[3][2]];
    [V_ci[4][1] V_ci[4][2]]]

plot_polygon(p_path, V_i, 2, :blue, "")
plot!(p_path, [s_hist[1][1]], [s_hist[1][2]], 
    markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="Initial State")

# plot final state
V_cf = pose_to_corners(s_hist[end], veh)
V_f = [[V_cf[1][1] V_cf[1][2]];
    [V_cf[2][1] V_cf[2][2]];
    [V_cf[3][1] V_cf[3][2]];
    [V_cf[4][1] V_cf[4][2]]]

plot_polygon(p_path, V_f, 2, :red, "")
plot!(p_path, [s_hist[end][1]], [s_hist[end][2]], 
    markercolor=:red, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="Final State")

display(p_path)