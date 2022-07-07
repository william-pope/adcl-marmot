# trajectory plot for physical system

using Plots
using BSON: @save, @load
using CSV
using DataFrames

include("../HJB-planner/HJB_generator_functions.jl")

hjb_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/HJB-planner/bson/"
hjb_path_nuc = "/home/adcl/Documents/marmot-algs/HJB-planner/bson/"
hjb_path = hjb_path_nuc

@load hjb_path*"env.bson" env
@load hjb_path*"veh.bson" veh

hist_path_mac = "/Users/willpope/Desktop/Research/marmot-ros/controller_pkg/histories/"
hist_path_nuc = "/home/adcl/catkin_ws/src/marmot-ros/controller_pkg/histories/"
hist_path = hist_path_nuc

# load planner histories
@load hist_path*"s_hist.bson" s_hist
veh_hist = []
ped1_hist = []
for s_k in s_hist
    push!(veh_hist, s_k[1:3])
    push!(ped1_hist, s_k[4:5])
end

# load VRPN histories
csv = CSV.File(hist_path*"veh_hist_vrpn.csv")
veh_hist_vrpn = []
for s_k in csv
    push!(veh_hist_vrpn, s_k)
end

csv = CSV.File(hist_path*"ped1_hist_vrpn.csv")
ped1_hist_vrpn = []
for s_kv in csv
    push!(ped1_hist_vrpn, s_kv)
end

@show size(veh_hist,1)
@show size(veh_hist_vrpn,1)

# plot optimal path from y_0 to target set
p_path = plot(aspect_ratio=:equal, size=(750,1050), 
            xlabel="x-axis [m]", ylabel="y-axis [m]",
            title="Vicon Output", 
            legend=false,
            titlefontsize = 20,
            legend_font_pointsize = 11,
            top_margin = -5*Plots.mm,
            left_margin = 4*Plots.mm)

plot_polygon(p_path, env.W, 2, :black, "Workspace")
plot_polygon(p_path, env.T_xy, 2, :green, "Target Set")
# plot_polygon(p_path, env.O_vec[1], 2, :red, "Obstacle")
# plot_polygon(p_path, env.O_vec[2], 2, :red, "")
# plot_polygon(p_path, env.O_vec[3], 2, :red, "")]

plot!(p_path, getindex.(veh_hist,1), getindex.(veh_hist,2).+7,
    linewidth=0, markershape=:circle, markersize=3, markerstrokewidth=0 ,
    label="Planning Points")

plot!(p_path, getindex.(veh_hist_vrpn,1), getindex.(veh_hist_vrpn,2).+7,
    linewidth=1,
    label="Vicon Output")

plot!(p_path, getindex.(ped1_hist,1), getindex.(ped1_hist,2).+7,
    linewidth=0, markershape=:circle, markersize=3, markerstrokewidth=0 ,
    label="Planning Points")

plot!(p_path, getindex.(ped1_hist_vrpn,1), getindex.(ped1_hist_vrpn,2).+7,
    linewidth=1,
    label="Vicon Output")

# plot initial state
V_ci = pose_to_edges(s_hist[1], veh)
V_i = [[V_ci[1][1] V_ci[1][2]];
    [V_ci[2][1] V_ci[2][2]];
    [V_ci[3][1] V_ci[3][2]];
    [V_ci[4][1] V_ci[4][2]]]

plot_polygon(p_path, V_i, 2, :blue, "")
plot!(p_path, [s_hist[1][1]], [s_hist[1][2]], 
    markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="Initial State")

# plot final state
V_cf = pose_to_edges(s_hist[end], veh)
V_f = [[V_cf[1][1] V_cf[1][2]];
    [V_cf[2][1] V_cf[2][2]];
    [V_cf[3][1] V_cf[3][2]];
    [V_cf[4][1] V_cf[4][2]]]

plot_polygon(p_path, V_f, 2, :red, "")
plot!(p_path, [s_hist[end][1]], [s_hist[end][2]], 
    markercolor=:red, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="Final State")

display(p_path)