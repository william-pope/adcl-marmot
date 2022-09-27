# HJB_plotting.jl

include("HJB_utils.jl")
include("dynamics_models.jl")

# TO-DO: need to fix to make work with 1-d value array
#   - may have to assemble n=dim value array, should be similar process to 
function plot_HJB_value(value_array, heatmap_clim, env, veh, sg)
    # reshape value array into n-dimensional array
    value_array_m = reshape(value_array, sg.state_grid.cut_counts...)

    # plot HJB value function as a heat map
    for i3_plot in eachindex(sg.state_grid.cutPoints[3])
        # for i4_plot in eachindex(sg.state_grid.cutPoints[4])
            val_plot = transpose(value_array_m[:, :, i3_plot])

            p_k = heatmap(sg.state_grid.cutPoints[1], sg.state_grid.cutPoints[2], 
                        val_plot, clim=(0, heatmap_clim),
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

            x = [x_pos, y_pos, sg.state_grid.cutPoints[3][i3_plot]]
            
            veh_body = state_to_body(x, veh)
                
            plot!(p_k, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

            # plot_polygon(p_k, V, 2, :blue, "Vehicle")
            plot!(p_k, veh_body, alpha=0.0, linecolor=:blue, linewidth=2, linealpha=1.0, label="Vehicle")

            plot!(p_k, [x_max], [y_pos], markercolor=:white, label="")
            plot!(p_k, [x_pos], [y_min], markercolor=:white, label="")

            theta_deg = round(rad2deg(x[3]), digits=1)
            annotate!(x_pos, y_pos+1.5, text("theta [deg]:\n$theta_deg", 14))

            display(p_k)
        # end
    end
end

# anim = @animate 
# gif(anim, algs_path*"HJB-planner/figures/hjb_theta.gif", fps=4)

function plot_HJB_growth(value_array, heatmap_clim, step, k_plot, env, veh)
    # reshape value array into n-dimensional array
    value_array_m = reshape(value_array, sg.state_grid.cut_counts...)
    
    # plot HJB value function as a heat map
    p_step = heatmap(sg.state_grid.cutPoints[1], sg.state_grid.cutPoints[2], 
                transpose(value_array_m[:,:,k_plot]), clim=(0, heatmap_clim),
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

    x = [x_pos, y_pos, sg.state_grid.cutPoints[3][k_plot]]
    
    veh_body = state_to_body(x, veh)
        
    plot!(p_step, [x_pos], [y_pos], markercolor=:blue, markershape=:circle, markersize=3, markerstrokewidth=0, label="")

    # plot_polygon(p_step, V, 2, :blue, "Vehicle")
    plot!(p_step, veh_body, alpha=0.0, linecolor=:blue, linewidth=2, linealpha=1.0, label="Vehicle")

    plot!(p_step, [x_max], [y_pos], markercolor=:white, label="")
    plot!(p_step, [x_pos], [y_min], markercolor=:white, label="")

    # plot step count
    # annotate!(x_pos, y_pos+1.5, text("step:\n$(step-1)", 14))

    display(p_step)
end

# anim_value_array = @animate 
# gif(anim_value_array, algs_path*"HJB-planner/figures/hjb_growth.gif", fps=3)

function plot_HJB_path(x_path_list)
    p_path = plot(aspect_ratio=:equal, 
                size=(800,800))

    plot!(p_path, env.workspace, alpha=0.0, linecolor=:black, linewidth=2, linealpha=1.0, label="Workspace")
    plot!(p_path, env.goal, alpha=0.0, linecolor=:green, linewidth=2, linealpha=1.0, label="Goal")

    if isempty(env.obstacle_list) == false
        plot!(p_path, env.obstacle_list[1], alpha=0.0, linecolor=:red, linewidth=2, linealpha=1.0, label="Obstacle")

        for obstacle in env.obstacle_list
            plot!(p_path, obstacle, alpha=0.0, linecolor=:red, linewidth=2, linealpha=1.0)
        end
    end

    for x_path in x_path_list
        # path
        plot!(p_path, getindex.(x_path,1), getindex.(x_path,2),
            # linez=getindex.(x_path,4),    
            linewidth = 2,
            label="")

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