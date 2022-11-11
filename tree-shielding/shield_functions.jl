# shield_functions.jl

# TO-DO: add Minkowski sum for human radius to each set
#   - Julia crashes when trying to use minkowski_sum() within shielding function
#   - gen_FRS() function works fine with mink_sum() when used separately in command line

# TO-DO: add finer time discretion using x_subpaths

function shield_action_set(x_k1, nearby_human_positions, Dt_obs_to_k1, Dt_plan, get_actions::Function, veh, human_goal_positions, v_human_max, m)
    Dv_max = 0.5    # NOTE: this should be pulled from one of the param structs

    # generate each human FRS sequence from t_k2 to t_stop_max
    actions_k1, ia_k1_set, _ = get_actions(x_k1, Dt_plan, veh, m)

    v_k2_max = x_k1[4] + maximum(getindex.(actions_k1, 2))
    kd_max = ceil(Int, (0.0 - v_k2_max)/(-Dv_max)) - 1

    F_all_body_seq = generate_F_all_seq(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_human_max, human_goal_positions, kd_max)

    # perform set check on all actions in standard POMDP action set
    ia_k1_safe_set = []

    for ia_k1 in ia_k1_set
        ia_k1_safe = false

        # propagate vehicle state to state x_k2
        a_k1 = actions_k1[ia_k1]
        x_k2, _ = propagate_state(x_k1, a_k1, Dt_plan, veh)

        # calculate time needed for divert path from new state
        kd_stop = ceil(Int, (0.0 - x_k2[4])/(-Dv_max)) - 1

        # iterate through divert steering angles
        for dpath in shuffle([1,2,3])
            dpath_safe = true

            x_kd = x_k2
            for kd in 0:kd_stop
                # TO-DO:
                # # check if divert path is in static environment RIC
                # val_x_kd = interp_value(x_kd, value_array, sg)
                # if val_x_kd < -50.0
                #   dpath_safe = false
                #   break
                # end

                # check for collisions with each human
                veh_body_cir_kd = state_to_body_circle(x_kd, veh)

                humans_safe = true
                for ih in axes(nearby_human_positions, 1)
                    F_ih_body_kd = F_all_body_seq[ih][3+kd]

                    if isdisjoint(veh_body_cir_kd, F_ih_body_kd) == false
                        humans_safe = false
                        break
                    end
                end

                if humans_safe == false
                    dpath_safe = false
                    break
                end

                # propagate vehicle to next step along divert path
                actions_kd, _, ia_divert_set = get_actions(x_kd, Dt_plan, veh, m)
                ia_d = ia_divert_set[dpath]
                a_d = actions_kd[ia_d]

                x_kd1, _ = propagate_state(x_kd, a_d, Dt_plan, veh)

                # pass state to next step
                x_kd = x_kd1
            end

            if dpath_safe == true
                ia_k1_safe = true
                break
            end
        end

        # action is safe
        if ia_k1_safe == true
            push!(ia_k1_safe_set, ia_k1)
        end
    end

    return actions_k1[ia_k1_safe_set], ia_k1_safe_set
end

# generate an FRS sequence for all nearby humans
function generate_F_all_seq(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_human_max, goal_positions, kd_max)
    F_all_body_seq = []

    # generate FRS sequences for each nearby human
    for x_ih_obs in nearby_human_positions
        F_ih_seq, F_ih_body_seq = generate_F_ih_seq(x_ih_obs, Dt_obs_to_k1, Dt_plan, v_human_max, goal_positions, kd_max)

        push!(F_all_body_seq, F_ih_body_seq)
    end

    return F_all_body_seq
end

# generate an FRS sequence for a given human position and time horizon
function generate_F_ih_seq(x_ih_obs, Dt_obs_to_k1, Dt_plan, v_human_max, goal_positions, kd_max)
    test = false

    F_ih_seq = []
    F_ih_body_seq = []

    h_body = VPolyCircle([0.0, 0.0], 0.381)

    # create initial polygon from observed position
    x_ih_ks_points = [x_ih_obs]
    F_ih_ks = VPolygon(x_ih_ks_points)
    push!(F_ih_seq, F_ih_ks)

    F_ih_body_ks = F_ih_ks
    # F_ih_body_ks = minkowski_sum(F_ih_ks, h_body)
    push!(F_ih_body_seq, F_ih_body_ks)

    # TEST ONLY ---
    if test == true
        p_F = plot(getindex.(goal_positions,1), getindex.(goal_positions,2), label="Goals",
            aspect_ratio=:equal, size=(800,800), dpi=300,
            linewidth=0,
            markershape=:circle, markersize=5,
            xticks=0:1:20, yticks=0:1:20)

        scatter!(p_F, getindex.(x_ih_ks_points, 1), getindex.(x_ih_ks_points, 2),
            markersize=2,
            label="")

        ks = 0
        annotate!(10.0+ks, 5.0, ks)

        display("image/png", p_F)

        plot!(p_F, F_ih_body_ks,
            label="")

        display("image/png", p_F)
    end
    # ---

    # propagate set through time steps
    for ks1 in 1:(2+kd_max)
        x_ih_ks1_points = Vector{Vector{Float64}}()

        ks1 == 1 ? Dt = Dt_obs_to_k1 : Dt = Dt_plan

        # apply all actions to each vertex of current set polygon
        for x_ih_ks in x_ih_ks_points
            for ig in axes(goal_positions, 1)
                x_ih_ks1 = propagate_human(x_ih_ks, ig, Dt, v_human_max, goal_positions)

                push!(x_ih_ks1_points, x_ih_ks1)
            end
        end

        # create set polygons from propagated states
        F_ih_ks1 = VPolygon(x_ih_ks1_points)
        push!(F_ih_seq, F_ih_ks1)

        F_ih_body_ks1 = minkowski_sum(F_ih_ks1, h_body)
        push!(F_ih_body_seq, F_ih_body_ks1)

        # TEST ONLY ---
        if test == true
            scatter!(p_F, getindex.(x_ih_ks1_points, 1), getindex.(x_ih_ks1_points, 2),
                markersize=2,
                label="")

            annotate!(10.0+ks1, 5.0, ks1)

            display("image/png", p_F)

            plot!(p_F, F_ih_body_ks1,
                label="")

            display("image/png", p_F)
        end
        # ---

        # pass states to next time step
        x_ih_ks_points = x_ih_ks1_points
    end

    return F_ih_seq, F_ih_body_seq
end

function propagate_human(x_ih_k, ig, Dt, v_human, goal_positions)
    # break out current state
    xp_ih_k = x_ih_k[1]
    yp_ih_k = x_ih_k[2]

    # pull out chosen goal location
    xpg = goal_positions[ig].x
    ypg = goal_positions[ig].y

    # calculate derivative at current state
    C_x = ((xpg-xp_ih_k)^2 + (ypg-yp_ih_k)^2)^(-1/2)

    xp_ih_dot_k = v_human * C_x * (xpg-xp_ih_k)
    yp_ih_dot_k = v_human * C_x * (ypg-yp_ih_k)

    # calculate next state
    xp_ih_k1 = xp_ih_k + (xp_ih_dot_k * Dt)
    yp_ih_k1 = yp_ih_k + (yp_ih_dot_k * Dt)

    # reassemble state vector
    x_ih_k1 = [xp_ih_k1, yp_ih_k1]

    return x_ih_k1
end

# local version for testing/plotting
function shield_action_set_test(x_k1, nearby_human_positions, Dt_obs_to_k1, Dt_plan, get_actions::Function, veh, human_goal_positions)
    test = false

    # TEST ONLY ---
    if test == true
        p1 = plot(aspect_ratio=:equal, size=(800,800), dpi=300,
            xticks=0:4:20, yticks=0:4:20,
            xlabel="x-axis [m]", ylabel="y-axis [m]",
            legend=:right)
    end
    # ---

    Dv_max = 0.5    # NOTE: this should be pulled from one of the param structs

    # generate each human FRS sequence from t_k2 to t_stop_max
    actions_k1, ia_k1_set, _ = get_actions(x_k1, Dt_plan, veh)

    v_k2_max = x_k1[4] + maximum(getindex.(actions_k1, 2))      # ISSUE: can reach over v_max=2.0 m/s limit
    kd_max = ceil(Int, (0.0 - v_k2_max)/(-Dv_max)) - 1

    v_human_max = 1.25
    F_all_body_seq = generate_F_all_seq(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_human_max, human_goal_positions, kd_max)

    # TEST ONLY ---
    if test == true
        println("x_k1 = ", x_k1)
        println("kd_max = ", kd_max)

        plot!(p1, [x_k1[1]], [x_k1[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        plot!(p1, state_to_body(x_k1, veh))
    end
    # ---

    # perform set check on all actions in standard POMDP action set
    ia_k1_safe_set = []

    for ia_k1 in ia_k1_set
        ia_k1_safe = false

        # propagate vehicle state to state x_k2
        a_k1 = actions_k1[ia_k1]
        x_k2, _ = propagate_state(x_k1, a_k1, Dt_plan, veh)     # ISSUE: can reach over v_max=2.0 m/s limit

        # calculate time needed for divert path from new state
        kd_stop = ceil(Int, (0.0 - x_k2[4])/(-Dv_max)) - 1

        # TEST ONLY ---
        if test == true
            println("\nia_k1 = ", ia_k1, ", a_k1 = ", a_k1)
            # println("x_k2 = ", x_k2)
            # println("kd_stop = ", kd_stop)
            plot!(p1, [x_k2[1]], [x_k2[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot!(p1, state_to_body(x_k2, veh))
        end
        # ---

        # iterate through divert steering angles
        for dpath in shuffle([1,2,3])
            dpath_safe = true

            # TEST ONLY ---
            if test == true
                println("dpath = ", dpath)
                x_path_divert = []
            end
            # ---

            x_kd = x_k2
            for kd in 0:kd_stop
                # TO-DO:
                # # check if divert path is in static environment RIC
                # val_x_kd = interp_value(x_kd, value_array, sg)
                # if val_x_kd < -50.0
                #   dpath_safe = false
                #   break
                # end

                # check for collisions with each human
                veh_body_cir_kd = state_to_body_circle(x_kd, veh)

                # TEST ONLY ---
                if test == true
                    println("kd = ", kd, ", x_k = ", x_kd)
                    push!(x_path_divert, x_kd)
                    plot!(p1, getindex.(x_path_divert, 1), getindex.(x_path_divert, 2), label="")
                    plot!(p1, veh_body_cir_kd)
                end
                # ---

                humans_safe = true
                for ih in axes(nearby_human_positions, 1)
                    F_ih_body_kd = F_all_body_seq[ih][3+kd]

                    # TEST ONLY ---
                    if test == true
                        plot!(p1, [nearby_human_positions[ih][1]], [nearby_human_positions[ih][2]], label="", markershape=:circle, markersize=5)
                        plot!(p1, F_ih_body_kd, label="")
                    end
                    # ---

                    if isdisjoint(veh_body_cir_kd, F_ih_body_kd) == false
                        humans_safe = false
                        break
                    end
                end

                # TEST ONLY ---
                if test == true
                    display(p1)
                end
                # ---

                if humans_safe == false
                    dpath_safe = false
                    break
                end

                # propagate vehicle to next step along divert path
                actions_kd, _, ia_divert_set = get_actions(x_kd, Dt_plan, veh)
                ia_d = ia_divert_set[dpath]
                a_d = actions_kd[ia_d]

                x_kd1, _ = propagate_state(x_kd, a_d, Dt_plan, veh)

                # TEST ONLY ---
                if test == true
                    println("a_d = ", a_d)
                end
                # ---

                # pass state to next step
                x_kd = x_kd1
            end

            if dpath_safe == true
                ia_k1_safe = true
                break
            end
        end

        # TEST ONLY ---
        if test == true
            println("ia_k1_safe = ", ia_k1_safe)
        end
        # ---

        # action is safe
        if ia_k1_safe == true
            push!(ia_k1_safe_set, ia_k1)
        end
    end

    return actions_k1[ia_k1_safe_set], ia_k1_safe_set
end