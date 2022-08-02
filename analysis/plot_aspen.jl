# trajectory plot for physical system

using Plots
using BSON: @save, @load
using CSV
using DataFrames

include("../HJB-planner/HJB_generator_functions.jl")
include("../HJB-planner/HJB_planner_functions.jl")

hjb_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/HJB-planner/bson/"
hjb_path_nuc = "/home/adcl/Documents/marmot-algs/HJB-planner/bson/"
hjb_path = hjb_path_mac

@load hjb_path*"env.bson" env
@load hjb_path*"veh.bson" veh

hist_path_mac = "/Users/willpope/Desktop/Research/marmot-ros/controller_pkg/histories/stored/"
hist_path_nuc = "/home/adcl/catkin_ws/src/marmot-ros/controller_pkg/histories/"
hist_path = hist_path_mac

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"
algs_path = algs_path_mac

# load planner histories
@load hist_path*"s_hist.bson" s_hist
veh_hist = []
# ped1_hist = []
for s_k in s_hist
    push!(veh_hist, s_k[1:3])
    # push!(ped1_hist, s_k[4:5])
end

# load VRPN histories
csv = CSV.File(hist_path*"veh_hist_vrpn.csv")
veh_hist_vrpn = []
for s_k in csv
    push!(veh_hist_vrpn, s_k)
end

# csv = CSV.File(hist_path*"ped1_hist_vrpn.csv")
# ped1_hist_vrpn = []
# for s_kv in csv
#     push!(ped1_hist_vrpn, s_kv)
# end

@show size(veh_hist,1)
@show size(veh_hist_vrpn,1)

@load algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
@load algs_path*"HJB-planner/bson/target_mat.bson" target_mat
@load algs_path*"HJB-planner/bson/obstacle_mat.bson" obstacle_mat
@load algs_path*"HJB-planner/bson/env.bson" env
@load algs_path*"HJB-planner/bson/veh.bson" veh

x_0 = veh_hist[1]
x_sim_path, u_path, step = HJB_planner(x_0, U_HJB, dt, plan_steps, A, obstacle_mat, env, veh)

# plot optimal path from y_0 to target set
p_path = plot(aspect_ratio=:equal, size=(605,1000), 
            # xlabel="x-axis [m]", ylabel="y-axis [m]",
            # title="Vicon Output", 
            legend=:topright,
            titlefontsize = 20,
            legend_font_pointsize = 11,
            top_margin = -5*Plots.mm,
            left_margin = 4*Plots.mm)

plot_polygon(p_path, env.W, 3, :black, "")
plot_polygon(p_path, env.T_xy, 3, :green, "")
plot_polygon(p_path, env.O_vec[1], 3, :red, "")
    for O in env.O_vec
        plot_polygon(p_path, O, 3, :red, "")
    end

# plot!(p_path, getindex.(veh_hist,1), getindex.(veh_hist,2),
#     linewidth=0, markershape=:circle, markersize=3, markerstrokewidth=0 ,
#     label="Planning Points")

plot!(p_path, getindex.(veh_hist_vrpn,1), getindex.(veh_hist_vrpn,2),
    linewidth=2, linecolor=:purple,
    label="Actual Path")

plot!(p_path, getindex.(x_sim_path,1), getindex.(x_sim_path,2),
    linewidth=2, linestyle=:dash,
    label="Simulated Path")

# plot!(p_path, getindex.(ped1_hist,1), getindex.(ped1_hist,2),
#     linewidth=0, markershape=:circle, markersize=3, markerstrokewidth=0 ,
#     label="Planning Points")

# plot!(p_path, getindex.(ped1_hist_vrpn,1), getindex.(ped1_hist_vrpn,2),
#     linewidth=1,
#     label="Vicon Output")

# plot initial state
V_ci = pose_to_edges(veh_hist_vrpn[1], veh)
V_i = [[V_ci[1][1] V_ci[1][2]];
    [V_ci[2][1] V_ci[2][2]];
    [V_ci[3][1] V_ci[3][2]];
    [V_ci[4][1] V_ci[4][2]]]

plot_polygon(p_path, V_i, 2, :blue, "")
plot!(p_path, [veh_hist_vrpn[1][1]], [veh_hist_vrpn[1][2]], 
    markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="")

# plot final state
V_cf = pose_to_edges(veh_hist_vrpn[end], veh)
V_f = [[V_cf[1][1] V_cf[1][2]];
    [V_cf[2][1] V_cf[2][2]];
    [V_cf[3][1] V_cf[3][2]];
    [V_cf[4][1] V_cf[4][2]]]

plot_polygon(p_path, V_f, 2, :red, "")
plot!(p_path, [veh_hist_vrpn[end][1]], [veh_hist_vrpn[end][2]], 
    markercolor=:red, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="")

# plot final state
V_cf = pose_to_edges(x_sim_path[end], veh)
V_f = [[V_cf[1][1] V_cf[1][2]];
    [V_cf[2][1] V_cf[2][2]];
    [V_cf[3][1] V_cf[3][2]];
    [V_cf[4][1] V_cf[4][2]]]

plot_polygon(p_path, V_f, 2, :red, "")
plot!(p_path, [x_sim_path[end][1]], [x_sim_path[end][2]], 
    markercolor=:red, markershape=:circle, markersize=3, markerstrokewidth=0, 
    label="")

display(p_path)