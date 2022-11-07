# unused.jl

function generate_human_FRS(x_human, v_human, Dt_obs_to_kd, goal_positions)
    reach_radius = v_human * Dt_obs_to_kd

    vertices = Vector{Vector{Float64}}(undef, length(goal_positions)+1)
    for ig in axes(goal_positions, 1)
        goal_vector = reach_radius * normalize(goal_positions[ig] - x_human)
        vertices[ig] = x_human + goal_vector
    end
    vertices[end] = x_human

    F_human = VPolygon(vertices)

    # TO-DO: add Minkowski sum for human radius
    # human_radius = 0.2
    # C_h = VPolyCircle([0.0, 0.0], human_radius)

    # @show typeof(F_human)
    # @show typeof(C_h)

    # F_hc = LazySets.minkowski_sum(F_human, C_h)

    # @show typeof(F_hc)
    
    return F_human
end

function generate_human_FRS_sequence(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_k2_max, Dv_max, v_human, goal_positions)
    # calculate steps for longest divert path
    kd_stop_max = ceil(Int, (0.0 - v_k2_max)/(-Dv_max))

    # println("\nkd_stop_max = ", kd_stop_max)

    # generate human reachable sets at each divert time step
    F_seq = []
    for kd = 0:kd_stop_max
        Dt_obs_to_kd = Dt_obs_to_k1 + Dt_plan + kd*Dt_plan
        # println("kd = ", kd, ", Dt_obs_to_kd = ", Dt_obs_to_kd)

        F_all_k = []
        for x_human in nearby_human_positions
            F_human_k = generate_human_FRS(x_human, v_human, Dt_obs_to_kd, goal_positions)

            push!(F_all_k, F_human_k)
        end

        push!(F_seq, F_all_k)
    end

    return F_seq
end