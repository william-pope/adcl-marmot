# Hamilton-Jacobi-Bellman demonstration

using Plots
using BSON: @save, @load
include("HJB_functions.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   
unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)   
veh = unit_car

# define workspace
# W = [[-3.010 -7.023];
#     [3.010 -7.023];
#     [3.010 7.023];
#     [-3.010 7.023]]

W = [[-4.0 -4.0];
    [4.0 -4.0];
    [4.0 4.0];
    [-4.0 4.0]]

# define target set
# T_xy = [[-0.5 2.0];
#         [0.5 2.0];
#         [0.5 3.5];
#         [-0.5 3.5]]

T_xy = [[-0.75 -0.75];
        [0.75 -0.75];
        [0.75 0.75];
        [-0.75 0.75]]

T_theta = [[-pi, pi]]

# define obstacles
O1 = [[-1.5 -2.5];
    [1.5 -2.5];
    [1.5 -2.0];
    [-1.5 -2.0]]

O2 = [[1.0 -1.0];
    [3.0 -0.5];
    [3.0 0.5];
    [1.0 0.0]]

O3 = [[-1.5 2.0];
    [1.5 2.0];
    [1.5 2.5];
    [-1.5 2.5]]

# O_vec = [O1]
O_vec = []

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
    max_reps = 5000
    anim_bool = true
    @time U_HJB = solve_HJB_PDE(du_tol, max_reps, env, veh, anim_bool)

    @save algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
    @save algs_path*"HJB-planner/bson/env.bson" env
    @save algs_path*"HJB-planner/bson/veh.bson" veh
else
    @load algs_path*"HJB-planner/bson/U_HJB.bson" U_HJB
end


# 4) PLOTS --- --- ---
veh = unit_car
x_max = 4 + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)
y_min = -5 - sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

# plot U as heat map
anim = @animate for k_plot in 1:size(env.theta_grid,1)
    p_k = heatmap(env.x_grid, env.y_grid, 
                transpose(U_HJB[:,:,k_plot]), clim=(0,10),
                aspect_ratio=:equal, size=(1050,1050),
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
    y = [2, -5, env.theta_grid[k_plot]]
    
    V_c = pose_to_corners(y, veh)
    V = [[V_c[1][1] V_c[1][2]];
        [V_c[2][1] V_c[2][2]];
        [V_c[3][1] V_c[3][2]];
        [V_c[4][1] V_c[4][2]]]
        
    plot!(p_k, [2], [y_min], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot!(p_k, [2], [-5], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot_polygon(p_k, V, 2, :blue, "Vehicle Orientation")

    theta_deg = round(rad2deg(y[3]), digits=1)
    annotate!(-2, -5, text("theta [deg]:\n$theta_deg", 14))

    # display(p_k)
end

gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)

