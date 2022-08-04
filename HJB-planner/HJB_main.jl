# Hamilton-Jacobi-Bellman demonstration

using Plots
using BSON: @save, @load
using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_planner_functions.jl")
include("HJB_utils.jl")
include("dynamics_models.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   
unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)   
veh = unit_car

EoM = dubins_car_EoM

# # Reeds-Shepp car
# actions = [[a_v,a_phi] for a_v in [-veh.u_vb_max, veh.u_vf_max], a_phi in [-veh.u_phi_max, 0.0, veh.u_phi_max]]

# Dubins car
actions = [[a_v,a_phi] for a_v in [veh.u_vf_max], a_phi in [-veh.u_phi_max, 0.0, veh.u_phi_max]]

actions = reshape(actions, (length(actions),1))
sort!(actions, dims=1)

# define workspace
# W = [[0.0 0.0];
#     [5.518 0.0];
#     [5.518 11.036];
#     [0.0 11.036]]

W = [[0.0 0.0];
    [8.0 0.0];
    [8.0 12.0];
    [0.0 12.0]]

# W = [[0.0 0.0];
#     [100.0 0.0];
#     [100.0 100.0];
#     [0.0 100.0]]

# define target set
# T_xy = [[3.15 9.8];
#         [4.35 9.8];
#         [4.35 11.0];
#         [3.15 11.0]]

T_xy = [[3.25 10.5];
        [4.75 10.5];
        [4.75 12.0];
        [3.25 12.0]]

# T_xy = [[95.0 83.0];
#         [100.0 83.0];
#         [100.0 87.0];
#         [95.0 87.0]]

T_theta = [[-pi, pi]]

# define obstacles
# - circular obstacles defined as [x, y, r], converted to polygon overapproximation
# OC1 = circle_to_polygon([1.5, 3.0, 0.6])
# OC2 = circle_to_polygon([4.0, 5.0, 0.6])
# OC3 = circle_to_polygon([2.8, 8.0, 0.6])

OC1 = circle_to_polygon([2.5, 3.0, 0.6])
OC2 = circle_to_polygon([6.0, 5.0, 0.6])
OC3 = circle_to_polygon([3.8, 9.0, 0.6])

# OC1 = circle_to_polygon([65.0, 35.0, 20.0])
# OC2 = circle_to_polygon([35.0, 70.0, 10.0])

O_vec = [OC1, OC2, OC3]
# O_vec = []

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

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"

algs_path = algs_path_mac

gen_HJB = false
plan_HJB = true
plot_result = true

# generate HJB value function
if gen_HJB == true
    du_tol = 0.01
    max_steps = 5000
    plot_growth = false
    # @btime solve_HJB_PDE(A, du_tol, max_steps, env, veh, anim_bool)
    value_array, target_array, obstacle_array = solve_HJB_PDE(actions, du_tol, max_steps, env, veh, EoM, plot_growth)

    N_grid = size(env.x_grid,1) * size(env.y_grid,1) * size(env.theta_grid,1)
    println("total grid nodes = ", N_grid)

    @save algs_path*"HJB-planner/bson/value_array.bson" value_array
    @save algs_path*"HJB-planner/bson/target_array.bson" target_array
    @save algs_path*"HJB-planner/bson/obstacle_array.bson" obstacle_array
    @save algs_path*"HJB-planner/bson/env.bson" env
    @save algs_path*"HJB-planner/bson/veh.bson" veh
else
    @load algs_path*"HJB-planner/bson/value_array.bson" value_array
    @load algs_path*"HJB-planner/bson/target_array.bson" target_array
    @load algs_path*"HJB-planner/bson/obstacle_array.bson" obstacle_array
    @load algs_path*"HJB-planner/bson/env.bson" env
    @load algs_path*"HJB-planner/bson/veh.bson" veh
end



# plan paths to goal
x_path_list = []

x_0_list = [[6.0, 1.5, deg2rad(90)],
            [6.75, 3.0, deg2rad(105)],
            [2.75, 5.2, deg2rad(-85)],
            [0.7, 8.0, deg2rad(25)]]

if plan_HJB == true
    dt = 0.2
    max_steps = 5e3

    for x_0 in x_0_list
        println("value at x_0 = ", interp_value_SL(x_0, value_array, env))

        x_path, u_path, step = plan_HJB_path(x_0, actions, dt, value_array, obstacle_array, max_steps, EoM, env, veh)
        push!(x_path_list, x_path)

        path_time = step*dt
        println("path execution time: ", path_time, " sec")
    end
end
    

# 4) PLOTS --- --- ---

if plot_result == true
    # plot U as heat map
    # anim = @animate 
    for k_plot in [28] #eachindex(env.theta_grid)
        p_k = heatmap(env.x_grid, env.y_grid, 
                    transpose(value_array[:,:,k_plot]), clim=(0,15),
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

        plot_polygon(p_k, env.W, 3, :black, "Workspace")
        plot_polygon(p_k, env.T_xy, 3, :green, "Target Set")
        plot_polygon(p_k, env.O_vec[1], 3, :red, "Obstacle")
        for O in env.O_vec
            plot_polygon(p_k, O, 3, :red, "")
        end

        # vehicle figure
        x_pos = 10.0
        y_pos = 4.5

        x_max = x_pos + sqrt((veh.ext_l-veh.ext2axle)^2 + (veh.ext_w/2)^2)
        y_min = y_pos - sqrt((veh.ext_l-veh.ext2axle)^2 + (veh.ext_w/2)^2)

        x = [x_pos, y_pos, env.theta_grid[k_plot]]
        
        V_c = pose_to_edges(x, veh)
        V = [[V_c[1][1] V_c[1][2]];
            [V_c[2][1] V_c[2][2]];
            [V_c[3][1] V_c[3][2]];
            [V_c[4][1] V_c[4][2]]]
            
        plot!(p_k, [x_pos], [y_pos], 
            markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

        plot_polygon(p_k, V, 2, :blue, "Vehicle")

        plot!(p_k, [x_max], [y_pos], markercolor=:white, label="")
        plot!(p_k, [x_pos], [y_min], markercolor=:white, label="")

        theta_deg = round(rad2deg(x[3]), digits=1)
        annotate!(x_pos, y_pos+1.5, text("theta [deg]:\n$theta_deg", 14))

        for x_path in x_path_list
            # path
            plot!(p_k, getindex.(x_path,1), getindex.(x_path,2),
                linewidth = 2, linecolor=:white,
                label="")
    
            # start position
            plot!(p_k, [x_path[1][1]], [x_path[1][2]], 
                markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, 
                label="")
    
            V_c = pose_to_edges(x_path[1], veh)
            V = [[V_c[1][1] V_c[1][2]];
                [V_c[2][1] V_c[2][2]];
                [V_c[3][1] V_c[3][2]];
                [V_c[4][1] V_c[4][2]]]
    
            plot_polygon(p_k, V, 2, :white, "")
    
            # end position
            plot!(p_k, [x_path[end][1]], [x_path[end][2]], 
                markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, 
                label="")
    
            V_c = pose_to_edges(x_path[end], veh)
            V = [[V_c[1][1] V_c[1][2]];
                [V_c[2][1] V_c[2][2]];
                [V_c[3][1] V_c[3][2]];
                [V_c[4][1] V_c[4][2]]]
    
            plot_polygon(p_k, V, 2, :white, "")
        end

        display(p_k)
    end

    # gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)
end

# @btime HJB_action(x_0, value_array, A, obstacle_array, env, veh)

# ProfileView.@profview for _ in 1:1000
#     HJB_action(x_0, value_array, A, obstacle_array, env, veh)
# end