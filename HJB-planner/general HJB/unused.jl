struct ActionGrid
    action_grid::RectangleGrid
    action_list_static::Array{Any}
end

# define action set
action_space = [[-0.5, 0.0, 0.5], [-0.475, 0.475]]
du_num_steps = [3, 5]
ag = define_action_grid(action_space, du_num_steps)

# discretizes action space
function define_action_grid(action_space, du_num_steps)
    action_iters = [range(minimum(axis), maximum(axis), du_num_steps[i]) for (i, axis) in enumerate(action_space)]
    action_grid = RectangleGrid(action_iters...)

    action_list_static = []
    for action in action_grid
        push!(action_list_static, SA[action...])
    end

    ag = ActionGrid(action_grid, action_list_static)
    return ag
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k, u_k, dt::Float64, EoM, veh, sg)    
    w1 = EoM(x_k, u_k, veh)
    w2 = EoM(x_k + w1*dt/2, u_k, veh)
    w3 = EoM(x_k + w2*dt/2, u_k, veh)
    w4 = EoM(x_k + w3*dt, u_k, veh)

    x_k1 = x_k + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    x_k1_m = map(mod_angle, x_k1, sg.angle_wrap_array) 

    return x_k1_m
end

function mod_angle(x::Float64, wrap::Bool)
    if wrap == true
        xm = x % (2*pi)
        if xm > pi
            xm -= 2*pi
        end
    else
        xm = x
    end

    return xm
end

# x_k1 = runge_kutta_4(x_k, a_k, Dt_plan, EoM, veh, sg)
# x_k1 = propagate_4d_v_state(x_k, a_k, Dt_plan, EoM, veh, sg)

# sortperm!(nbr_indices, nbr_weights, rev=true) 
# @show nbr_indices

# sort!(nbr_indices, by = nbr_indices[x] -> nbr_weights[x])
# @show nbr_indices

# Dv = value_k1 - value_k
# println(ia)
# println("s->a->sp (iter): ", x_k, " -> ", ag.action_grid[ia], " -> ", x_k1)
# println("v->vp: ", value_k, " -> ", value_k1, ", Dv = ", Dv)

# println("taking threshold action: ", ia, "\n")
# println("taking fully minimized action: ", ia_min, "\n")

# TO-DO: improve action selection method
#   - goal: check as few actions as possible, while still returning a reasonable action
#   - for all possible actions, could create some ordered list of their values
#       - to be fast, want to choose action near top of list, without actually calculating list
#   - know:
#       - optimal action at each neighboring node
#       - change in value for an optimal action (equal to Dt_plan)
#   - near-optimal value cutoff:
#       - in order to avoid checking whole action space for best value, need some cutoff condition to end search
#       - cutoff threshold is an estimate, could have rare case that none of the actions meet your condition
#           - in this case, will just have to search full action space and take best available
#   - ordering: 
#       - want to order indices of surrounding nodes by their distance to the state
#       - remaining action indices filled in at end, randomized
#       - no idea how to do this lol, needs to be efficient somehow

# STATUS: seems like it works...
# TO-DO: test for Dubin's car, more x_0 points

# ISSUE: struggles to get exactly into goal region
#   - think that epsilon might be too forgiving, causing it to take ugly actions
#   - when flow field converges near goal, neighboring nodes more likely to have different actions
#   - (?): does the optimal action always have Dv=-Dt_plan? think interpolation/solving method get in way of this
#   - fixes:
#       - increasing epsilon (more restrictive) makes it better, but still see loops for some paths
#       - just increasing the goal works lol
#       - could tighten epsilon as you get closer to goal (use straightline dist or value knowledge)

# x_k = SA[4.1, 1.7, deg2rad(120)]

# @btime HJB_policy(x_k, Dt, value_array, EoM, veh, sg)
# @btime fast_policy(x_k, Dt, value_array, opt_ia_array, EoM, veh, sg)

# ProfileView.@profview for _ in 1:1000
#     fast_policy(x_k, Dt, value_array, opt_ia_array, EoM, veh, sg)
# end

# X) PERFORMANCE --- --- ---


# x_dot = [0.0, 0.0, 0.0]
# x = [3.1, 1.2, 0.7*pi]
# a = [1.0, -0.475]

# x_dot = MVector(0.0, 0.0, 0.0)
# x = SVector(3.1, 1.23, 0.7*pi)
# a = SVector(1.0, -0.475)

# ProfileView.@profview 

# @btime bicycle_3d_EoM(x, a, veh)
# @btime bicycle_3d_EoM!(x_dot, x, a, veh)

# println("\nEoM")
# @btime bicycle_3d_EoM($x, $a, $veh)

# println("\nRK4")
# @btime runge_kutta_4($x, $a, $Dt, $EoM, $veh, $sg)

# println("\ninterp")
# @btime interpolate($sg.state_grid, $value_array, $x)
# @btime interp_value($x, $value_array, $sg)

# println("\npolicy (", length(actions), " actions)")
# @btime HJB_policy($x, $Dt, $value_array, $EoM, $veh, $sg)
# @btime fast_policy($x, $Dt, $value_array, $opt_ia_array, $EoM, $veh, $sg)

# @profview 
# fast_policy(x, Dt, value_array, opt_ia_array, EoM, veh, sg)

# Dt = 0.5
# substeps = 4

# x_k = SA[2.0, 3.5, deg2rad(0), 1.5]
# a_k = [0.0, 0.327]

# println("x_k = ", x_k)
# println("a_k = ", a_k)

# x_k1s = common_4d_DT_EoM_HJB(x_k, a_k, Dt)

# println("\nsingle step")
# println("x_k1 = ", x_k1)

# x_k1, x_k1_subpath = common_prop_HJB(x_k, a_k, Dt, substeps)

# println("\nsubpath")
# println("x_k1 = ", x_k1)
# println("x_k1_subpath = ", x_k1_subpath)

# p1 = plot(aspect_ratio=:equal)

# plot!(p1, [x_k[1]], [x_k[2]], markershape=:circle, markersize=5)
# plot!(p1, state_to_body(x_k, veh))

# plot!(p1, [x_k1s[1]], [x_k1s[2]], markershape=:circle, markersize=5)
# plot!(p1, state_to_body(x_k1s, veh))

# plot!(p1, getindex.(x_k1_subpath,1), getindex.(x_k1_subpath,2), markershape=:circle, markersize=5)
# plot!(p1, state_to_body(x_k1, veh))

# display(p1)