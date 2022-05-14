# Hamilton-Jacobi-Bellman demonstration

using Plots
using BSON: @save, @load
include("HJB_functions.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   
# unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)   
veh = marmot

# define workspace
W = [[-3.010 -7.023];
        [3.010 -7.023];
        [3.010 7.023];
        [-3.010 7.023]]

# define target set
T_xy = [[-0.5 -0.5];
            [0.5 -0.5];
            [0.5 0.5];
            [-0.5 0.5]]

T_theta = [[-pi, pi]]

# define obstacles
O1 = [[-3.0 -2.0];
        [-1.5 -2.5];
        [-0.5 -2.5];
        [-1.0 -1.0];
        [-3.0 -1.0]]

O2 = [[1.0 -1.0];
        [3.0 -0.5];
        [3.0 0.5];
        [1.0 0.0]]

O3 = [[-1.5 2.0];
        [1.5 2.0];
        [1.5 2.5];
        [-1.5 2.5]]

# O_vec = [O1, O2, O3]
O_vec = []

# initialize state grid
h_xy = 0.1
h_theta = deg2rad(5)

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

# ISSUE: value estimate not matching actual path length (time)
#   - is this ok/expected? look into further

# calculate HJB
run_HJB = false

if run_HJB == true
    du_tol = 0.01
    max_reps = 5000
    @time U_HJB = solve_HJB_PDE(du_tol, max_reps, env, veh, false)

    @save "/home/adcl/Documents/marmot-algs/HJB-planner/bson/U_HJB.bson" U_HJB
    @save "/home/adcl/Documents/marmot-algs/HJB-planner/bson/env.bson" env
    @save "/home/adcl/Documents/marmot-algs/HJB-planner/bson/veh.bson" veh
else
    @load "/home/adcl/Documents/marmot-algs/HJB-planner/bson/U_HJB.bson" U_HJB
end

# plan path
dt_path = 0.05


# 4) PLOTS --- --- ---
veh = marmot
x_max = 4 + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

# plot U as heat map
anim = @animate for k_plot in 1:size(env.theta_grid,1)
    p_k = heatmap(env.x_grid, env.y_grid, 
                transpose(U_HJB[:,:,k_plot]), clim=(0,10),
                aspect_ratio=:equal, size=(750,1050),
                xlabel="x-axis [m]", ylabel="y-axis [m]", 
                title="HJB Value Function",
                titlefontsize = 20,
                colorbar_title = "time-to-target [s]",
                legend=:topright,
                legend_font_pointsize = 11,
                top_margin = -30*Plots.mm,
                left_margin = 8*Plots.mm,
                bottom_margin = -8*Plots.mm)

    plot_polygon(p_k, env.W, 2, :black, "Workspace")
    plot_polygon(p_k, env.T_xy, 2, :green, "Target Set")
    # plot_polygon(p_k, env.O_vec[1], 2, :red, "Obstacle")
    # plot_polygon(p_k, env.O_vec[2], 2, :red, "")
    # plot_polygon(p_k, env.O_vec[3], 2, :red, "")

    # vehicle figure
    y = [4, -4, env.theta_grid[k_plot]]
    
    V_c = pose_to_corners(y, veh)
    V = [[V_c[1][1] V_c[1][2]];
        [V_c[2][1] V_c[2][2]];
        [V_c[3][1] V_c[3][2]];
        [V_c[4][1] V_c[4][2]]]
        
    plot!(p_k, [x_max], [-4], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot!(p_k, [4], [-4], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot_polygon(p_k, V, 2, :blue, "Vehicle Orientation")

    theta_deg = round(rad2deg(y[3]), digits=1)
    annotate!(4.25, 0, text("theta [deg]:\n$theta_deg", 14))

    display(p_k)
end

# gif(anim, "hjb_theta.gif", fps=4)

aspen_poses = [[-2, -6.5, 3*pi/4],
            [-1, -6.5, -pi/4],
            [0, -5.5, pi/2],
            [1, -6.5, -3*pi/4],
            [2, -6.5, pi/4]]

ex_poses = [[-4, -4, 0],
                [1, -4, pi/2],
                [3.5, -3.5, -pi/2]]

initial_poses = [[0.0, -5.5, pi/2]]

# plot optimal path from y_0 to target set
p_path = plot(aspect_ratio=:equal, size=(750,1050), 
            xlabel="x-axis [m]", ylabel="y-axis [m]",
            title="HJB Path Planner", 
            titlefontsize = 20,
            legend_font_pointsize = 11,
            top_margin = -5*Plots.mm,
            left_margin = 4*Plots.mm)

plot_polygon(p_path, env.W, 2, :black, "Workspace")
plot_polygon(p_path, env.T_xy, 2, :green, "Target Set")
# plot_polygon(p_path, env.O_vec[1], 2, :red, "Obstacle")
# plot_polygon(p_path, env.O_vec[2], 2, :red, "")
# plot_polygon(p_path, env.O_vec[3], 2, :red, "")

for y_0 in initial_poses
    (y_path, u_path) = HJB_planner(y_0, U_HJB, dt_path, env, veh)
    plot!(p_path, getindex.(y_path,1), getindex.(y_path,2), linewidth=2, label="")

    V_c = pose_to_corners(y_0, veh)
    V = [[V_c[1][1] V_c[1][2]];
        [V_c[2][1] V_c[2][2]];
        [V_c[3][1] V_c[3][2]];
        [V_c[4][1] V_c[4][2]]]

    plot!(p_path, [y_0[1]], [y_0[2]], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot_polygon(p_path, V, 2, :blue, "")
end

display(p_path)