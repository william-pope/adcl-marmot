
using Plots
include("HJB_functions.jl")

# 1) DEFINITIONS --- --- ---

function lidar_scan(s_veh, b_env, s_env_true, O_set, H_set, veh::Vehicle)
    radial_max = 10.0
    dr = 1.0
    angular_max = pi/2
    dpsi = deg2rad(10)

    b_env_p = deepcopy(b_env)
    
    lidar_base = [s_veh[1] + veh.wb/2*cos(s_veh[3]),
                s_veh[2] + veh.wb/2*sin(s_veh[3])]
    
    p_scan = plot(aspect_ratio=:equal)

    for psi in -angular_max : dpsi : angular_max
        for r in 0 : dr : radial_max
            pt = [lidar_base[1] + r*cos(s_veh[3] + psi),
                lidar_base[2] + r*sin(s_veh[3] + psi)]

            next_psi = false

            # checks if point is in a known obstacle
            for Oi_set in O_set
                if point_in_polygon(pt, Oi_set) == true
                    next_psi = true
                    break
                end
            end

            # checks if point is in spot for a hidden obstacle
            for i in 1:size(H_set,1)
                Hi_set = H_set[i]
                if point_in_polygon(pt, Hi_set) == true
                    # updates belief based on if spot is/isn't occupied
                    if s_env_true[i] == 1
                        b_env_p[i] = 1
                    else
                        b_env_p[i] = 0
                    end
                
                    next_psi = true
                    break
                end
            end

            if next_psi == true
                break
            end

            plot!(p_scan, [pt[1]], [pt[2]], markershape=:circle, markersize=4, label="")
        end
    end

    display(p_scan)

    return b_env_p
end

function point_in_polygon(pt, P_set)
    P_rows = [collect(1:size(P_set, 1)); 1]

    P_ineqs = zeros(Bool, size(P_set, 1))
    for i in 1:size(P_set, 1)
        i1 = P_rows[i]
        i2 = P_rows[i+1]

        x1 = P_set[i1,1]
        y1 = P_set[i1,2]
        x2 = P_set[i2,1]
        y2 = P_set[i2,2]

        val = (y1 - y2)*pt[1] + (x2 - x1)*pt[2] + x1*y2 - x2*y1

        if val >= 0
            P_ineqs[i] = 1
        end
    end

    if all(P_ineqs) == true
        return true
    else
        return false
    end
end


# 2) PARAMETERS --- --- ---

# vehicle parameters 
truck = Vehicle(5.0, 2.5, 0.5, 4.0, 6.0, 2.0, 1.0)  # pickup truck

# define workspace
W_set = [[-50.0 -50.0];
        [50.0 -50.0];
        [50.0 50.0];
        [-50.0 50.0]]

W_x_bounds = [minimum(W_set[:,1]), maximum(W_set[:,1])]
W_y_bounds = [minimum(W_set[:,2]), maximum(W_set[:,2])]

# define initial pose
y_0 = [-48.75, -44.0, 0]

# define target set
T_xy_set = [[40.0 42.5];
            [50.0 42.5];
            [50.0 47.5];
            [40.0 47.5]]

T_theta_set = [[deg2rad(-180), deg2rad(-175)],
                [deg2rad(175), deg2rad(180)]]

# define known obstacles
O1_set = [[-50.0 -40.0];
        [-35.0 -40.0];
        [-35.0 -38.0];
        [-50.0 -38.0]]

O2_set = [[-50.0 -50.0];
        [-35.0 -50.0];
        [-35.0 -48.0];
        [-50.0 -48.0]]

O3_set = [[35.0 38.0];
        [50.0 38.0];
        [50.0 40.0];
        [35.0 40.0]]

O4_set = [[35.0 48.0];
        [50.0 48.0];
        [50.0 50.0];
        [35.0 50.0]]

O_set = [O1_set, O2_set, O3_set, O4_set]

# define hidden obstacle locations and occupancy
H1_set = [[-5.0 10.0];
        [5.0 10.0];
        [5.0 12.5];
        [-5.0 12.5]]

H2_set = [[5.0 0.0];
        [10.0 0.0];
        [10.0 3.0];
        [5.0 3.0]]

H_set = [H1_set, H2_set]

s_env_true = [1, 1]

# initialize state grid
h_xy = 2.0
h_theta = deg2rad(8)

sg = StateGrid(h_xy, 
                h_theta,
                minimum(W_set[:,1]) : h_xy : maximum(W_set[:,1]),
                minimum(W_set[:,2]) : h_xy : maximum(W_set[:,2]),
                -pi : h_theta : pi)


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

s_veh = [0, 0, 0]
b_env = [0.5, 0.5]

b_env_p = lidar_scan(s_veh, b_env, s_env_true, O_set, H_set)

# du_tol = 0.01
# max_reps = 100
# @time U = solve_HJB_PDE(du_tol, max_reps, sg, O_set, truck)

# dt = 0.1

# (y_path, u_path) = find_path(y_0, U, T_xy_set, T_theta_set, dt, sg, truck)

# # ISSUE: still having problem with chattering (grid size might be only solution)

# # 4) PLOTS --- --- ---

# # plot U as heat map
# for k_plot in LinRange(1, size(sg.theta_grid, 1), 5)
#     k_plot = Int(round(k_plot, digits=0))
#     theta_k = round(rad2deg(sg.theta_grid[k_plot]), digits=3)

#     p_k = heatmap(sg.x_grid, sg.y_grid, transpose(U[:,:,k_plot]), clim=(0,100),
#                 aspect_ratio=:equal, size=(700,700), legend=:bottomright,
#                 xlabel="x-axis", ylabel="y-axis", title="HJB Value Function: u(x, y, theta=$theta_k)")

#     plot_polygon(W_set, p_k, 2, :black, "Workspace")
#     plot_polygon(T_xy_set, p_k, 2, :green, "Target Set")
#     plot_polygon(O1_set, p_k, 2, :red, "Obstacle")
#     plot_polygon(O2_set, p_k, 2, :red, "")

#     display(p_k)
# end


# # plot optimal path from y_0 to target set
# p_path = plot(aspect_ratio=:equal, size=(700,700), legend=:bottomright,
#             title="HJB Path Planner", xlabel="x-axis [m]", ylabel="y-axis [m]")

# plot_polygon(W_set, p_path, 2, :black, "Workspace")
# plot_polygon(T_xy_set, p_path, 2, :green, "Target Set")
# plot_polygon(O1_set, p_path, 2, :red, "Obstacle")
# plot_polygon(O2_set, p_path, 2, :red, "")

# plot!(p_path, getindex.(y_path,1), getindex.(y_path,2), linewidth=2, label="")

# display(p_path)