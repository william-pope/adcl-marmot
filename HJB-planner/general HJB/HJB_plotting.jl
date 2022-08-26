# HJB_plotting.jl

include("HJB_utils.jl")
include("dynamics_models.jl")

function plot_HJB_result(value_array, x_path_list, env, veh)
    # plot HJB value function as a heat map
    for k_plot in eachindex(sg.grid_array[3])
        p_k = heatmap(sg.grid_array[1], sg.grid_array[2], 
                    transpose(value_array[:,:,k_plot]), clim=(0,150),
                    aspect_ratio=:equal, 
                    size=(800,800),
                    # xlabel="x-axis [m]", ylabel="y-axis [m]", 
                    # title="HJB Value Function",
                    titlefontsize = 20,
                    legend=false, 
                    # legend=:topright,
                    legend_font_pointsize = 11,
                    colorbar=false,
                    colorbar_title = "time-to-target [s]",
                    top_margin = -30*Plots.mm,
                    bottom_margin = 4*Plots.mm,
                    left_margin = 8*Plots.mm,
                    right_margin = 0*Plots.mm)

        # TO-DO: plot as outlines with no fill
        plot!(p_k, env.workspace, color="blue", alpha=0.25, label="Workspace")
        plot!(p_k, env.goal, color="green", alpha=0.25, label="Goal")

        # if isempty(env.O_vec) == false
        #     plot_polygon(p_k, env.O_vec[1], 3, :red, "Obstacle")
        #     for O in env.O_vec
        #         plot_polygon(p_k, O, 3, :red, "")
        #     end
        # end

        # vehicle figure
        x_pos = sg.grid_array[1][end] + 1.5
        y_pos = sg.grid_array[2][end]/2  - 0.5

        x_max = x_pos + sqrt((veh.axis_to_cent_x + 1/2*veh.body_length)^2 + (veh.axis_to_cent_y + 1/2*veh.body_width)^2)
        y_min = y_pos - sqrt((veh.axis_to_cent_x + 1/2*veh.body_length)^2 + (veh.axis_to_cent_y + 1/2*veh.body_width)^2)

        x = [x_pos, y_pos, sg.grid_array[3][k_plot]]
        
        veh_body = state_to_body(x, veh)
            
        plot!(p_k, [x_pos], [y_pos], 
            markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

        # plot_polygon(p_k, V, 2, :blue, "Vehicle")
        plot!(p_k, veh_body, color="blue", label="Vehicle")

        plot!(p_k, [x_max], [y_pos], markercolor=:white, label="")
        plot!(p_k, [x_pos], [y_min], markercolor=:white, label="")

        # theta_deg = round(rad2deg(x[3]), digits=1)
        # annotate!(x_pos, y_pos+1.5, text("theta [deg]:\n$theta_deg", 14))

        # # plot paths over heat map
        # for x_path in x_path_list
        #     # path
        #     plot!(p_k, getindex.(x_path,1), getindex.(x_path,2),
        #         linewidth = 2, linecolor=:white,
        #         label="")

        #     # start position
        #     plot!(p_k, [x_path[1][1]], [x_path[1][2]], 
        #         markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, 
        #         label="")

        #     V_c = pose_to_edges(x_path[1], veh)
        #     V = [[V_c[1][1] V_c[1][2]];
        #         [V_c[2][1] V_c[2][2]];
        #         [V_c[3][1] V_c[3][2]];
        #         [V_c[4][1] V_c[4][2]]]

        #     plot_polygon(p_k, V, 2, :white, "")

        #     # end position
        #     plot!(p_k, [x_path[end][1]], [x_path[end][2]], 
        #         markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, 
        #         label="")

        #     V_c = pose_to_edges(x_path[end], veh)
        #     V = [[V_c[1][1] V_c[1][2]];
        #         [V_c[2][1] V_c[2][2]];
        #         [V_c[3][1] V_c[3][2]];
        #         [V_c[4][1] V_c[4][2]]]

        #     plot_polygon(p_k, V, 2, :white, "")
        # end

        display(p_k)
    end
end

# anim = @animate 
# gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)

function plot_HJB_growth(value_array, step, theta_plot, env::Environment, veh::Vehicle)
    k_plot = find_idx(theta_plot, env.theta_grid) + 1

    p_step = heatmap(env.x_grid, env.y_grid, transpose(value_array[:,:,k_plot]), clim=(0,15),
                # xlim=(-3.5,5.5),
                aspect_ratio=:equal, 
                size=(750,1000),
                # xlabel="x-axis [m]", ylabel="y-axis [m]", 
                colorbar_title = "time-to-target [s]",
                legend=:topright,
                # legend=true, 
                colorbar=false,
                legend_font_pointsize = 11,
                top_margin = -30*Plots.mm,
                bottom_margin = 4*Plots.mm,
                left_margin = 8*Plots.mm,
                right_margin = 0*Plots.mm)

    plot_polygon(p_step, env.W, 3, :black, "Workspace")
    plot_polygon(p_step, env.T_xy, 3, :green, "Target Set")
    if isempty(env.O_vec) == false
        plot_polygon(p_step, env.O_vec[1], 3, :red, "Obstacle")
        for O in env.O_vec
            plot_polygon(p_step, O, 3, :red, "")
        end
    end

    # plot vehicle figure
    x_pos = env.x_grid[end] + 1.0
    y_pos = env.y_grid[end]/2

    x_max = x_pos + sqrt((veh.axle_l-veh.ext2axle)^2 + (veh.ext_w/2)^2)
    y_min = y_pos - sqrt((veh.axle_l-veh.ext2axle)^2 + (veh.ext_w/2)^2)

    x = [x_pos, y_pos, env.theta_grid[k_plot]]
    
    E_arr = pose_to_edges(x, veh)
    V = [[E_arr[1][1] E_arr[1][2]];
        [E_arr[2][1] E_arr[2][2]];
        [E_arr[3][1] E_arr[3][2]];
        [E_arr[4][1] E_arr[4][2]]]
        
    plot!(p_step, [x_max], [y_pos], markercolor=:white, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot!(p_step, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot_polygon(p_step, V, 2, :blue, "Vehicle")

    # plot step count
    annotate!(x_pos, y_pos+1.5, text("step:\n$(step-1)", 14))

    display(p_step)
end

# anim_value_array = @animate 
# gif(anim_value_array, algs_path*"HJB-planner/figures/hjb_growth.gif", fps=3)