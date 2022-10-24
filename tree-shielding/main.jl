# main.jl

# overall objective:
#   - given a vehicle state and pedestrian positions, return the actions that are guaranteed safe

# need to collision-check for:
#   - every next state sp (~10)
#       - every divert trajectory (3)
#           - every time step (~5)
#               - every pedestrian set (~6)

# assumed parameters:
#   - vehicle top speed: 3 m/s
#   - vehicle acceleration: +/- 1.5 m/s
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

        for k in time_steps
            propagate vehicle state to (t+k*dt)

            for human in near_humans
                generate F_p(t+k*dt)
                check polygon intersection between vehicle and human set
            end

            check HJB value to see if in static obstacle or RIC
        end
    end
end
==#

# F_cir = Ball2(human, v_h*t_hrz)
# plot!(p1, F_cir, alpha=0.0, linecolor=:blue, linestyle=:dash, linewidth=1, linealpha=1)

using LazySets
using LinearAlgebra
using Random
using Plots
using BenchmarkTools

include("utils.jl")

# define environment
env_width = 5.5
env_length = 11.0

goals = [[0.0, 0.0],
        [env_width, 0.0],
        [env_width, env_length],
        [0.0, env_length]]


# define vehicle and dynamics
wheelbase = 0.324
body_dims = [0.5207, 0.2762]
origin_to_cent = [0.1715, 0.0]
veh = define_vehicle(wheelbase, body_dims, origin_to_cent)

EoM = bicycle_4d_v_EoM


# x_n1 = [2.0, 1.0, pi/2, 2.0]

# plot environment and state x_0
p1 = plot(getindex.(goals,1), getindex.(goals,2), label="Goals",
    aspect_ratio=:equal, size=(500,600), linewidth=0, 
    markershape=:circle, markersize=5)



# plot!(p1, [x_n1[1]], [x_n1[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
# plot!(p1, state_to_body(x_n1, veh))


function generate_human_FRS(x_h, v_h, t_horizon, goals)
    reach_radius = v_h * t_horizon

    vertices = Vector{Vector{Float64}}(undef, length(goals)+1)
    for i in eachindex(goals)
        goal_vector = reach_radius * normalize(goals[i] - x_h)
        vertices[i] = x_h + goal_vector
    end
    vertices[end] = x_h

    F_h = VPolygon(vertices)

    human_radius = 0.2
    C_h = VPolyCircle([0.0, 0.0], human_radius)

    # @show typeof(F_h)
    # @show typeof(C_h)

    # F_hc = LazySets.minkowski_sum(F_h, C_h)

    # @show typeof(F_hc)
    
    return F_h
end


# TO-DO: see if more efficient on human polygons to check all paths at same time step at once
#   - might lose some of the gains in the sequential checking logic

# TO-DO: check that k_divert is accurate, especially for velocities/time steps that don't line up nicely

# TO-DO: check HJB value at each step in divert path

# ISSUE: Julia crashes when trying to use minkowski_sum() within shielding function
#   - gen_FRS() function works fine with mink_sum() when used separately in command line

function shield_action_set(x_0, actions, Dt_plan, Dt_divert, v_max, EoM, veh)
    # main shielding function
    max_brake = -2.0
    divert_angles = [-phi_max, 0.0, phi_max]

    # calculate maximum time horizon
    t_horizon = (0 - v_max) / max_brake
    k_horizon = round(Int, t_horizon / Dt_divert)

    # generate human reachable sets at each divert time step
    F_hist = []
    for k = 0:k_horizon
        t_k = k*Dt_divert

        F_k_list = []
        for human in humans
            F_human_k = generate_human_FRS(human, v_h, 2*Dt_plan+t_k, goals)

            # return 1
            push!(F_k_list, F_human_k)
        end

        push!(F_hist, F_k_list)
    end

    safe_actions = []

    num_checks = 0

    plot!(p1, [x_0[1]], [x_0[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
    plot!(p1, state_to_body(x_0, veh))

    # iterate through all actions in standard POMDP action set at current time step
    for ia in eachindex(actions)
        # propagate vehicle state to next time step
        x_1 = runge_kutta_4(x_0, actions[ia], Dt_plan, EoM, veh)   

        x_1_valid = false

        # calculate time needed for divert path from new state
        t_divert = (0 - x_1[4]) / max_brake
        k_divert = round(Int, t_divert / Dt_divert)

        println("\na_0 = ", actions[ia])
        println("x_1 = ", x_1)
        println("t_divert = ", t_divert)
        println("k_divert = ", k_divert)

        plot!(p1, [x_1[1]], [x_1[2]], markercolor=:black, markershape=:circle, markersize=3, markerstrokewidth=0, label="")
        plot!(p1, state_to_body(x_1, veh))

        # iterate through divert steering angles
        for phi_d in shuffle(divert_angles)
            a_d = [max_brake, phi_d]

            divert_path = []
            divert_valid = true
            
            x_k = x_1
            for k in 0:k_divert
                # propagate vehicle along divert path to stop
                t_k = k*Dt_divert
                veh_body_k = state_to_body(x_k, veh)

                println("k = ", k, ", x_k = ", x_k)

                push!(divert_path, x_k)
                plot!(p1, getindex.(divert_path, 1), getindex.(divert_path, 2), label="")
                plot!(p1, veh_body_k)

                # TO-DO: check HJB value of x_k for ststic obstacle/RIC
                
                # iterate through each human
                for h in 1:length(humans)
                    num_checks += 1

                    # collision-check
                    if isempty(intersection(veh_body_k, F_hist[k+1][h])) == false
                        # println("collision at k = ", k, ", x_k = ", x_k)
                        divert_valid = false
                    end

                    plot!(p1, [humans[h][1]], [humans[h][2]], label="", markershape=:circle, markersize=5)
                    plot!(p1, F_hist[k+1][h], label="")

                    # step in collision
                    if divert_valid == false
                        break
                    end
                end

                display(p1)

                # path in collision
                if divert_valid == false
                    break

                # path is safe
                elseif divert_valid == true && k == k_divert
                    x_1_valid = true

                end

                x_k1 = runge_kutta_4(x_k, a_d, Dt_divert, EoM, veh) 
                x_k = x_k1
            end

            # action is safe
            if x_1_valid == true
                push!(safe_actions, ia)
                break
            end
        end
    end

    println("num_checks = ", num_checks)
    return safe_actions
end

# define human positions and velocity
humans = [[2.0, 7.5],
        [3.8, 5.2],
        [4.3, 8.4]]
        # ,
        # [1.6, 3.2],
        # [2.8, 5.1],
        # [4.8, 4.6]]

v_h = 1.0

# define action set
v_max = 2.5
phi_max = 0.35
actions = [[0.0, -phi_max],
            [0.0, 0.0],
            [0.0, phi_max]]
            
            # ,   
            # [0.0, -1/3*phi_max],  
            # [0.0, 1/3*phi_max], 
            # [0.0, -2/3*phi_max],  
            # [0.0, 2/3*phi_max],  
            # [-1.0, 0.0],
            # [1.0, 0.0]]

Dt_plan = 0.5
Dt_divert = 0.25

x_0 = [2.5, 1.2, pi/2, 2.0]
        
safe_actions = shield_action_set(x_0, actions, Dt_plan, Dt_divert, v_max, EoM, veh)

# @btime shield_action_set($x_0, $actions, $Dt_plan, $Dt_divert, $v_max, $EoM, $veh)

# performance:
#   - for 9 actions, 6 humans, no optimization
#       - 9/9 safe: 310.040 us, 4792 allocations
#       - 4/9 safe: 463.547 us, 7070 allocations
#       - 0/9 safe: 192.423 us, 3069 allocations