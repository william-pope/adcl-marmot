# HJB_plotting.jl

using Plots

include("HJB_utils.jl")

# plot final value array over all velocity and heading states
function plot_HJB_value(value_array, heatmap_clim, env, veh, sg)
    # reshape value array into n-dimensional array
    value_array_m = reshape(value_array, sg.state_grid.cut_counts...)

    # plot HJB value function as a heat map
    for i3_plot in eachindex(sg.state_grid.cutPoints[3])
        for i4_plot in eachindex(sg.state_grid.cutPoints[4])
            plot_val = transpose(value_array_m[:, :, i3_plot, i4_plot])

            p_k = heatmap(sg.state_grid.cutPoints[1], sg.state_grid.cutPoints[2], 
                        plot_val, clim=(0, heatmap_clim),
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

            plot!(p_k, env.workspace, alpha=0.0, linecolor=:black, linewidth=2, linealpha=1.0, label="Workspace")
            plot!(p_k, env.goal, alpha=0.0, linecolor=:green, linewidth=2, linealpha=1.0, label="Goal")

            if isempty(env.obstacle_list) == false
                plot!(p_k, env.obstacle_list[1], alpha=0.0, linecolor=:red, linewidth=2, linealpha=1.0, label="Obstacle")

                for obstacle in env.obstacle_list
                    plot!(p_k, obstacle, alpha=0.0, linecolor=:red, linewidth=2, linealpha=1.0)
                end
            end

            # vehicle figure
            x_pos = sg.state_grid.cutPoints[1][end] + 1.5
            y_pos = sg.state_grid.cutPoints[2][end]/2  - 0.5

            x_max = x_pos + sqrt((veh.origin_to_cent[1] + 1/2*veh.body_dims[1])^2 + (veh.origin_to_cent[2] + 1/2*veh.body_dims[2])^2)
            y_min = y_pos - sqrt((veh.origin_to_cent[1] + 1/2*veh.body_dims[1])^2 + (veh.origin_to_cent[2] + 1/2*veh.body_dims[2])^2)

            x = [x_pos, y_pos, sg.state_grid.cutPoints[3][i3_plot], sg.state_grid.cutPoints[4][i4_plot]]
            
            veh_body = state_to_body(x, veh)
                
            plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

            plot!(p_k, veh_body, alpha=0.0, linecolor=:blue, linewidth=2, linealpha=1.0, label="Vehicle")

            plot!(p_k, [x_max], [y_pos], markercolor=:white, label="")
            plot!(p_k, [x_pos], [y_min], markercolor=:white, label="")

            theta_deg = round(rad2deg(x[3]), digits=1)
            annotate!(x_pos, y_pos+1.5, text("theta [deg]:\n$theta_deg", 14))

            v = round(x[4], digits=2)
            annotate!(x_pos, y_pos-1.5, text("velocity [m/s]:\n$v", 14))

            display(p_k)
        end
    end
end

# plot current value array at a given velocity/heading at each step in the solving process
function plot_HJB_growth(value_array, heatmap_clim, step, env, veh)
    # reshape value array into n-dimensional array
    value_array_m = reshape(value_array, sg.state_grid.cut_counts...)

    i3_plot = 19
    i4_plot = 1
    plot_val = transpose(value_array_m[:, :, i3_plot, i4_plot])
    
    # plot HJB value function as a heat map
    p_step = heatmap(sg.state_grid.cutPoints[1], sg.state_grid.cutPoints[2], 
                plot_val, clim=(0, heatmap_clim),
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

    plot!(p_step, env.workspace, alpha=0.0, linecolor=:black, linewidth=2, linealpha=1.0, label="Workspace")
    plot!(p_step, env.goal, alpha=0.0, linecolor=:green, linewidth=2, linealpha=1.0, label="Goal")

    if isempty(env.obstacle_list) == false
        plot!(p_step, env.obstacle_list[1], alpha=0.0, linecolor=:red, linewidth=2, linealpha=1.0, label="Obstacle")

        for obstacle in env.obstacle_list
            plot!(p_step, obstacle, alpha=0.0, linecolor=:red, linewidth=2, linealpha=1.0)
        end
    end

    # vehicle figure
    x_pos = sg.state_grid.cutPoints[1][end] + 1.5
    y_pos = sg.state_grid.cutPoints[2][end]/2  - 0.5

    x_max = x_pos + sqrt((veh.origin_to_cent[1] + 1/2*veh.body_dims[1])^2 + (veh.origin_to_cent[2] + 1/2*veh.body_dims[2])^2)
    y_min = y_pos - sqrt((veh.origin_to_cent[1] + 1/2*veh.body_dims[1])^2 + (veh.origin_to_cent[2] + 1/2*veh.body_dims[2])^2)

    x = [x_pos, y_pos, sg.state_grid.cutPoints[3][i3_plot], sg.state_grid.cutPoints[4][i4_plot]]
    
    veh_body = state_to_body(x, veh)
        
    plot!(p_step, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

    plot!(p_step, veh_body, alpha=0.0, linecolor=:blue, linewidth=2, linealpha=1.0, label="Vehicle")

    plot!(p_step, [x_max], [y_pos], markercolor=:white, label="")
    plot!(p_step, [x_pos], [y_min], markercolor=:white, label="")

    # step count
    annotate!(x_pos, y_pos+1.5, text("step:\n$(step-1)", 14))

    display(p_step)
end

# plot paths from planner
function plot_HJB_path(x_path_list, x_subpath_list)
    p_path = plot(aspect_ratio=:equal, 
                size=(800,800),
                xlabel="X-axis [m]", ylabel="Y-axis [m]",
                title="Rollout Policy vs Pure HJB Policy")

    # workspace
    plot!(p_path, env.workspace, 
        alpha=0.0, 
        linecolor=:black, linewidth=2, linealpha=1.0, 
        label="Workspace")
    
    # goal
    plot!(p_path, env.goal, 
        color=:green, alpha=0.125, 
        linecolor=:green, linewidth=2, linealpha=1.0, 
        label="Goal")

    # obstacles
    if isempty(env.obstacle_list) == false
        plot!(p_path, env.obstacle_list[1], 
            color=:red, alpha=0.125, 
            linecolor=:red, linewidth=2, linealpha=1.0, 
            label="Obstacle")

        for obstacle in env.obstacle_list[2:end]
            plot!(p_path, obstacle, 
                color=:red, alpha=0.125, 
                linecolor=:red, linewidth=2, linealpha=1.0,
                label="")
        end
    end

    label_list = ["Rollout Policy", "Pure HJB Policy"]
    for ip in 1:length(x_path_list)
        x_path = x_path_list[ip]
        x_subpath = x_subpath_list[ip]

        # shift velocity up one step to make line_z look right
        linez_velocity = zeros(length(x_subpath))
        for kk in 1:(length(x_subpath)-1)
            linez_velocity[kk] = x_subpath[kk+1][4]
        end

        # subpath lines
        plot!(p_path, getindex.(x_subpath,1), getindex.(x_subpath,2),
            linez=linez_velocity, clim=(0,3.5), colorbar_title="Velocity [m/s]",
            linewidth = 2,
            # markershape=:circle, markersize=1.5, markerstrokewidth=0, 
            label="")

        # path points
        plot!(p_path, getindex.(x_path,1), getindex.(x_path,2),
        linewidth = 0, linealpha=0.0,
        markershape=:circle, markersize=2.5, markerstrokewidth=0, 
        label=label_list[ip])

        # start position
        plot!(p_path, [x_path[1][1]], [x_path[1][2]], 
            markershape=:circle, markersize=3, markerstrokewidth=0, 
            label="")

        veh_body = state_to_body(x_path[1], veh)
        plot!(p_path, veh_body, alpha=0.0,
            linewidth=2, linealpha=1.0, label="")

        # end position
        plot!(p_path, [x_path[end][1]], [x_path[end][2]], 
            markershape=:circle, markersize=3, markerstrokewidth=0, 
            label="")

        veh_body = state_to_body(x_path[end], veh)
        plot!(p_path, veh_body, alpha=0.0, 
            linewidth=2, linealpha=1.0, label="")
    end

    display(p_path)
end

# anim = @animate 
# gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)

# anim_value_array = @animate 
# gif(anim_value_array, algs_path*"HJB-planner/figures/hjb_growth.gif", fps=3)