# main.jl

# overall notes:
#   - objective: given a vehicle state and pedestrian positions, return the actions that are guaranteed safe
#   - formulate for after-DESPOT approach
#   - don't want to change POMDP action set (for now)
#       - divert paths should be pulling from the same action parameters
#   - can approximate t_obs ~= t_k1 to simplify implementation in simulator

# integration:
#   - inputs: human observation, Dt_obs_to_k1, vehicle state

# need to collision-check for:
#   - every next state sp (~10)
#       - every divert trajectory (3)
#           - every time step (~5)
#               - every pedestrian set (~6)

# assumed parameters:
#   - vehicle top speed: 2.0 m/s (this shouldn't matter)
#   - vehicle max brake: -0.5 m/s (same as POMDP) (although would make sense to include as a harder brake?)
#   - human speed: 1 m/s

# pseudo-code:
#==
function gen_ped_bound(s0, t_h)
    return VPolygon
end

NOTE: better order to do this

for action in action_set
    propagate vehicle forward one time step from s to sp

    for path in divert_paths
        action_k = [full brake, steering]

        for k in time_steps (prefer dt=0.1 sec)
            propagate vehicle state to (t+k*dt)

            for human in near_humans
                generate F_p(t+k*dt)
                check polygon intersection between vehicle and human set (for current time step)
            end

            check HJB value to see if in static obstacle or RIC (better to convert HJB RICs to polygons over-approx to check definitively?)
        end
    end
end
==#

#=
notional architecture:

while end_run == false
    publish action -> a_k

    retrieve state/observation -> obs_k
    update belief -> belief_k

    ranked actions = run_DESPOT(belief_k, predicted_x_k1)
    [... 0.4 sec ...]

    retrieve another state/observation -> obs_shield
    take current time -> Dt_obs_to_k1
    ia_k1_safe_set = run_shield(obs_shield, Dt_obs_to_k1)

    a_k1 = best safe action
    a_k = a_k1
end

(?): can DESPOT.jl return Q-value for all actions?
=#

# performance:
#   - for 9 actions, 6 humans, no optimization
#       - 9/9 safe: 310.040 us, 4792 allocationss
#       - 4/9 safe: 463.547 us, 7070 allocations
#       - 0/9 safe: 192.423 us, 3069 allocations

#   - updated version -> still <500 us


# TO-DO: need to fix Minkwoski sum issue
# F_cir = Ball2(human, v_human*t_hrz)
# plot!(p1, F_cir, alpha=0.0, linecolor=:blue, linestyle=:dash, linewidth=1, linealpha=1)

using LazySets
using LinearAlgebra
using StaticArrays
using Random
using Plots
using BenchmarkTools

using Pkg
Pkg.develop(PackageSpec(path = "/Users/willpope/.julia/dev/BellmanPDEs"))
using BellmanPDEs

# define environment
ws_width = 5.5
ws_length = 11.0

goal_positions = [[0.0, 0.0],
                [ws_width, 0.0],
                [ws_width, ws_length],
                [0.0, ws_length]]

# define vehicle and dynamics
wheelbase = 0.324
body_dims = [0.5207, 0.2762]
origin_to_cent = [0.1715, 0.0]
phi_max = 0.475
v_max = 2.0
discount = 1.0
veh = define_vehicle(wheelbase, body_dims, origin_to_cent, phi_max, v_max, discount)

# define action set
function get_actions(x, Dt, veh)
    # set change in velocity (Dv) limit
    Dv_lim = 0.5

    # set steering angle (phi) limit
    Dtheta_lim = deg2rad(45)

    v = x[4]
    vp = v + Dv_lim
    vn = v - Dv_lim

    phi_lim = atan(Dtheta_lim * 1/Dt * 1/abs(v) * veh.l)
    phi_lim = clamp(phi_lim, 0.0, veh.phi_max)

    phi_lim_p = atan(Dtheta_lim * 1/Dt * 1/abs(vp) * veh.l)
    phi_lim_p = clamp(phi_lim_p, 0.0, veh.phi_max)

    phi_lim_n = atan(Dtheta_lim * 1/Dt * 1/abs(vn) * veh.l)
    phi_lim_n = clamp(phi_lim_n, 0.0, veh.phi_max)

    actions = SVector{9, SVector{2, Float64}}(
        (-phi_lim_n, -Dv_lim),        # Dv = -Dv
        (0.0, -Dv_lim),
        (phi_lim_n, -Dv_lim),

        (-phi_lim, 0.0),       # Dv = 0.0
        (0.0, 0.0),
        (phi_lim, 0.0),

        (-phi_lim_p, Dv_lim),        # Dv = +Dv
        (0.0, Dv_lim),
        (phi_lim_p, Dv_lim))

    ia_set = collect(1:length(actions))

    return actions, ia_set
end

# plot environment and state x_0
p1 = plot(getindex.(goal_positions,1), getindex.(goal_positions,2), label="Goals",
    aspect_ratio=:equal, size=(500,600), linewidth=0, 
    markershape=:circle, markersize=5,
    xticks=0:1:6, yticks=0:1:10)



# TO-DO: modify to use more general set calculation
#   - shouldn't affect larger architecture
#   - may need to pass in set from previous time step in order to propagate
#   - easier to generate sequence for one human at a time (similar to numerically propagating vehicle forward along subpath)

function propagate_human(x_human_k, ig, Dt, v_human, goal_positions)
    # break out current state
    xph_k = x_human_k[1]
    yph_k = x_human_k[2]

    # pull out chosen goal location
    xpg = goal_positions[ig][1]
    ypg = goal_positions[ig][2]

    # calculate derivative at current state
    C_x = ((xpg-xph_k)^2 + (ypg-yph_k)^2)^(-1/2)

    xph_dot_k = v_human * C_x * (xpg-xph_k)
    yph_dot_k = v_human * C_x * (ypg-yph_k)

    # calculate next state
    xph_k1 = xph_k + (xph_dot_k * Dt)
    yph_k1 = yph_k + (yph_dot_k * Dt)

    # reassemble state vector
    x_human_k1 = [xph_k1, yph_k1]

    return x_human_k1
end

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


# ISSUE: Julia crashes when trying to use minkowski_sum() within shielding function
#   - gen_FRS() function works fine with mink_sum() when used separately in command line


function shield_action_set(x_k, ia_k, nearby_human_positions, Dt_obs_to_k1, get_actions::Function, Dt_plan, veh)
    # define divert actions (hard-coded)
    ia_divert_set = [1, 2, 3]
    Dv_max = 0.5    # NOTE: this should be pulled from one of the param structs
    
    # propagate vehicle to state x_k1
    actions_k, _ = get_actions(x_k, Dt_plan, veh)
    
    a_k = actions_k[ia_k]
    x_k1, _ = propagate_state(x_k, a_k, Dt_plan, veh)

    # println("x_k1 = ", x_k1)

    # generate human FRS sequence from t_k2 to t_stop_max
    actions_k1, ia_k1_set = get_actions(x_k1, Dt_plan, veh)
    v_k2_max = x_k1[4] + maximum(getindex.(actions_k1, 2))

    F_seq = generate_human_FRS_sequence(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_k2_max, Dv_max, v_human, goal_positions)

    # TEST ONLY ---
    plot!(p1, [x_k1[1]], [x_k1[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot!(p1, state_to_body(x_k1, veh))
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
        # println("\na_k1 = ", a_k1)
        # println("x_k2 = ", x_k2)
        # println("kd_stop = ", kd_stop)
        plot!(p1, [x_k2[1]], [x_k2[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        plot!(p1, state_to_body(x_k2, veh))
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
                # println("kd = ", kd, ", x_k = ", x_kd)
                push!(divert_path, x_kd)
                plot!(p1, getindex.(divert_path, 1), getindex.(divert_path, 2), label="")
                plot!(p1, veh_body_kd)
                # ---
                
                humans_safe = true
                for ih in axes(nearby_human_positions, 1)
                    
                    # TEST ONLY ---
                    plot!(p1, [nearby_human_positions[ih][1]], [nearby_human_positions[ih][2]], label="", markershape=:circle, markersize=5)
                    plot!(p1, F_seq[kd+1][ih], label="")
                    # ---

                    if isempty(intersection(veh_body_kd, F_seq[kd+1][ih])) == false || isempty(intersection(F_seq[kd+1][ih], veh_body_kd)) == false
                        # println("collision at kd = ", kd, ", x_kd = ", x_kd)
                        humans_safe = false
                        break
                    end
                end

                # TEST ONLY ---
                display(p1)
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


# define human positions and velocity
nearby_human_positions = [[2.0, 7.5],
                        [3.8, 5.2],
                        [4.3, 8.4],
                        [4.6, 3.2],
                        [2.8, 5.1],
                        [4.8, 4.6]]

# (?): what upper bound velocity should be used here?
v_human = 1.0

Dt_plan = 0.5

x_k = SVector(2.0, 1.0, pi*2/3, 1.5)
ia_k = 5
Dt_obs_to_k1 = 0.0
   
# ia_k1_safe_set = shield_action_set(x_k, ia_k, nearby_human_positions, Dt_obs_to_k1, get_actions, Dt_plan, veh)

# @btime shield_action_set($x_k, $ia_k, $nearby_human_positions, $Dt_obs_to_k1, $get_actions, $Dt_plan, $Dv_max, $veh)


x_human_k = nearby_human_positions[2]
ig = 1
Dt = 0.5

x_human_k1 = propagate_human(x_human_k, ig, Dt, v_human, goal_positions)