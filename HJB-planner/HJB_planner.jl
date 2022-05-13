# Hamilton-Jacobi-Bellman demonstration

using Plots
using JLD
include("HJB_functions.jl")

# 2) PARAMETERS --- --- ---

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   
# unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)   

# define workspace
W_set = [[-3.010 -7.023];
        [3.010 -7.023];
        [3.010 7.023];
        [-3.010 7.023]]

# define target set
T_xy_set = [[-0.5 4.0];
            [0.5 4.0];
            [0.5 5.0];
            [-0.5 5.0]]

T_theta_set = [[-pi, pi]]

# define obstacles
O1_set = [[-3.0 -2.0];
        [-1.5 -2.5];
        [-0.5 -2.5];
        [-1.0 -1.0];
        [-3.0 -1.0]]

O2_set = [[1.0 -1.0];
        [3.0 -0.5];
        [3.0 0.5];
        [1.0 0.0]]

O3_set = [[-1.5 2.0];
        [1.5 2.0];
        [1.5 2.5];
        [-1.5 2.5]]

# O_set = [O1_set, O2_set, O3_set]
O_set = []

# initialize state grid
h_xy = 0.1
h_theta = deg2rad(5)

sg = StateGrid(h_xy, 
                h_theta,
                minimum(W_set[:,1]) : h_xy : maximum(W_set[:,1]),
                minimum(W_set[:,2]) : h_xy : maximum(W_set[:,2]),
                -pi : h_theta : pi)


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

# ISSUE: value estimate not matching actual path length (time)
#   - is this ok/expected? look into further

# calculate HJB
run_HJB = false

if run_HJB == true
    du_tol = 0.01
    max_reps = 5000
    @time U_HJB = solve_HJB_PDE(du_tol, max_reps, sg, O_set, marmot, false)

    save("U_HJB.jld", "U_HJB", U_HJB)
    save("sg.jld", "sg", sg)
    save("veh.jld", "veh", marmot)
else
    U_HJB = load("U_HJB.jld", "U_HJB")
end

# plan path
dt_path = 0.05


# 4) PLOTS --- --- ---
veh = marmot
x_max = 4 + sqrt((veh.l-veh.b2a)^2 + (veh.w/2)^2)

# plot U as heat map
anim = @animate for k_plot in 1:size(sg.theta_grid,1)
    p_k = heatmap(sg.x_grid, sg.y_grid, 
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

    plot_polygon(W_set, p_k, 2, :black, "Workspace")
    plot_polygon(T_xy_set, p_k, 2, :green, "Target Set")
    # plot_polygon(O1_set, p_k, 2, :red, "Obstacle")
    # plot_polygon(O2_set, p_k, 2, :red, "")
    # plot_polygon(O3_set, p_k, 2, :red, "")

    # vehicle figure
    y = [4, -4, sg.theta_grid[k_plot]]
    
    veh_set = pose_to_corners(y, veh)
    veh_mat = [[veh_set[1][1] veh_set[1][2]];
                [veh_set[2][1] veh_set[2][2]];
                [veh_set[3][1] veh_set[3][2]];
                [veh_set[4][1] veh_set[4][2]]]
        
    plot!(p_k, [x_max], [-4], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot!(p_k, [4], [-4], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot_polygon(veh_mat, p_k, 2, :blue, "Vehicle Orientation")

    theta_deg = round(rad2deg(y[3]), digits=1)
    annotate!(4.25, 0, text("theta [deg]:\n$theta_deg", 14))

    display(p_k)
end

# gif(anim, "hjb_theta.gif", fps=4)

aspen_poses = [[-2, -6.5, 3*pi/4],
            [-1, -6.5, -pi/4],
            [0, -6.5, pi/2],
            [1, -6.5, -3*pi/4],
            [2, -6.5, pi/4]]

ex_poses = [[-4, -4, 0],
                [1, -4, pi/2],
                [3.5, -3.5, -pi/2]]

initial_poses = aspen_poses

# plot optimal path from y_0 to target set
p_path = plot(aspect_ratio=:equal, size=(750,1050), 
            xlabel="x-axis [m]", ylabel="y-axis [m]",
            title="HJB Path Planner", 
                titlefontsize = 20,
                legend_font_pointsize = 11,
            top_margin = -5*Plots.mm,
                left_margin = 4*Plots.mm)

plot_polygon(W_set, p_path, 2, :black, "Workspace")
plot_polygon(T_xy_set, p_path, 2, :green, "Target Set")
# plot_polygon(O1_set, p_path, 2, :red, "Obstacle")
# plot_polygon(O2_set, p_path, 2, :red, "")
# plot_polygon(O3_set, p_path, 2, :red, "")

for y_0 in initial_poses
    (y_path, u_path) = HJB_planner(y_0, U_HJB, T_xy_set, T_theta_set, dt_path, sg, veh)
    plot!(p_path, getindex.(y_path,1), getindex.(y_path,2), linewidth=2, label="")

    veh_set = pose_to_corners(y_0, veh)
    veh_mat = [[veh_set[1][1] veh_set[1][2]];
                [veh_set[2][1] veh_set[2][2]];
                [veh_set[3][1] veh_set[3][2]];
                [veh_set[4][1] veh_set[4][2]]]

    plot!(p_path, [y_0[1]], [y_0[2]], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot_polygon(veh_mat, p_path, 2, :blue, "")
end

display(p_path)