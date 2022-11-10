# shield_functions.jl

# TO-DO: add Minkowski sum for human radius to each set
#   - Julia crashes when trying to use minkowski_sum() within shielding function
#   - gen_FRS() function works fine with mink_sum() when used separately in command line

# TO-DO: get correct v_human to use (should be upper bound)

# TO-DO: add finer time discretion using x_subpaths

# TO-DO: configure inputs for POMDP integration
#   - create separate copy of code

function shield_action_set(x_k, ia_k, nearby_human_positions, Dt_obs_to_k1, Dt_plan, get_actions::Function, veh)
    test = true

    # TEST ONLY ---
    if test == true
        p1 = plot(getindex.(goal_positions,1), getindex.(goal_positions,2), label="Goals",
            aspect_ratio=:equal, size=(500,600), linewidth=0, 
            markershape=:circle, markersize=5,
            xticks=0:1:6, yticks=0:1:10)
            
        plot!(p1, [x_k[1]], [x_k[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        plot!(p1, state_to_body(x_k, veh))
    end
    # ---
    
    # define divert actions (hard-coded)
    ia_divert_set = [1, 2, 3]
    Dv_max = 0.5    # NOTE: this should be pulled from one of the param structs
    
    # propagate vehicle to state x_k1
    actions_k, _ = get_actions(x_k, Dt_plan, veh)
    
    a_k = actions_k[ia_k]
    x_k1, _ = propagate_state(x_k, a_k, Dt_plan, veh)

    # generate human FRS sequence from t_k2 to t_stop_max
    actions_k1, ia_k1_set = get_actions(x_k1, Dt_plan, veh)
    v_k2_max = x_k1[4] + maximum(getindex.(actions_k1, 2))
    kd_max = ceil(Int, (0.0 - v_k2_max)/(-Dv_max))

    # println("kd_max = ", kd_max)

    F_all_body_seq = generate_F_all_seq(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_human, goal_positions, kd_max)

    # TEST ONLY ---
    if test == true
        # println("x_k1 = ", x_k1)

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
        x_k2, _ = propagate_state(x_k1, a_k1, Dt_plan, veh)

        # calculate time needed for divert path from new state
        kd_stop = ceil(Int, (0.0 - x_k2[4])/(-Dv_max))

        # TEST ONLY ---
        if test == true
            # println("\na_k1 = ", a_k1)
            # println("x_k2 = ", x_k2)
            # println("kd_stop = ", kd_stop)
            plot!(p1, [x_k2[1]], [x_k2[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
            plot!(p1, state_to_body(x_k2, veh))
        end
        # ---

        # iterate through divert steering angles
        for ia_d in shuffle(ia_divert_set)
            ia_d_safe = true

            # TEST ONLY ---
            divert_path = []
            # ---

            x_kd = x_k2
            for kd in 0:kd_stop
                # TO-DO:
                # # check if divert path is in static environment RIC
                # val_x_kd = interp_value(x_kd, value_array, sg)
                # if val_x_kd < -50.0
                #   ia_d_safe = false
                #   break
                # end

                # check for collisions with each human
                veh_body_kd = state_to_body(x_kd, veh)

                # TEST ONLY ---
                if test == true
                    # println("kd = ", kd, ", x_k = ", x_kd)
                    push!(divert_path, x_kd)
                    plot!(p1, getindex.(divert_path, 1), getindex.(divert_path, 2), label="")
                    plot!(p1, veh_body_kd)
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

                    # if isempty(intersection(veh_body_kd, F_ih_body_kd)) == false || isempty(intersection(F_ih_body_kd, veh_body_kd)) == false
                    #     # println("collision at kd = ", kd, ", x_kd = ", x_kd)
                    #     humans_safe = false
                    #     break
                    # end

                    if isdisjoint(veh_body_kd, F_ih_body_kd) == false
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
                    ia_d_safe = false
                    break
                end
                
                # propagate vehicle to next step along divert path
                actions_kd, _ = get_actions(x_kd, Dt_plan, veh)
                a_d = actions_kd[ia_d]

                x_kd1, _ = propagate_state(x_kd, a_d, Dt_plan, veh)

                # println("a_d = ", a_d)

                # pass state to next step
                x_kd = x_kd1
            end

            if ia_d_safe == true
                ia_k1_safe = true
                break
            end 
        end

        # println("ia_k1_safe = ", ia_k1_safe)

        # action is safe
        if ia_k1_safe == true
            push!(ia_k1_safe_set, ia_k1)
        end
    end

    return ia_k1_safe_set
end

function generate_F_all_seq(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_human, goal_positions, kd_max)
    F_all_body_seq = []

    # generate FRS sequences for each nearby human
    for x_ih_obs in nearby_human_positions
        F_ih_seq, F_ih_body_seq = generate_F_ih_seq(x_ih_obs, Dt_obs_to_k1, Dt_plan, v_human, goal_positions, kd_max)
        
        push!(F_all_body_seq, F_ih_body_seq)
    end

    return F_all_body_seq
end

function generate_F_ih_seq(x_ih_obs, Dt_obs_to_k1, Dt_plan, v_human, goal_positions, kd_max)
    F_ih_seq = []
    F_ih_body_seq = []
    
    h_body = VPolyCircle([0.0, 0.0], 0.381)

    # create initial polygon from observed position
    x_ih_ks_points = [x_ih_obs]
    F_ih_ks = VPolygon(x_ih_ks_points)
    push!(F_ih_seq, F_ih_ks)

    F_ih_body_ks = minkowski_sum(F_ih_ks, h_body)
    push!(F_ih_body_seq, F_ih_body_ks)

    # propagate set through time steps
    for ks1 in 1:(2+kd_max)
        x_ih_ks1_points = Vector{Vector{Float64}}()

        ks1 == 1 ? Dt = Dt_obs_to_k1 : Dt = Dt_plan

        # apply all actions to each vertex of current set polygon 
        for x_ih_ks in x_ih_ks_points
            for ig in axes(goal_positions, 1)
                x_ih_ks1 = propagate_human(x_ih_ks, ig, Dt, v_human, goal_positions)
                
                push!(x_ih_ks1_points, x_ih_ks1)
            end
        end

        # create set polygons from propagated states
        F_ih_ks1 = VPolygon(x_ih_ks1_points)
        push!(F_ih_seq, F_ih_ks1)

        F_ih_body_ks1 = minkowski_sum(F_ih_ks1, h_body)
        push!(F_ih_body_seq, F_ih_body_ks1)

        # pass states to next time step
        x_ih_ks_points = x_ih_ks1_points
    end

    return F_ih_seq, F_ih_body_seq
end