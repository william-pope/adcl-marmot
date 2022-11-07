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
#   - every action a_k1 (~10)
#       - every divert trajectory (3)
#           - every time step (~5)
#               - every pedestrian set (~6)

# assumed parameters:
#   - vehicle max brake: -0.5 m/s (same as POMDP) (although would make sense to include a harder brake?)
#   - human speed: 1.0 m/s

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

#   - with generalized set prop -> 2.15 ms, 10499 allocations

using LazySets
using LinearAlgebra
using StaticArrays
using Random
using Plots
using BenchmarkTools

using Pkg
Pkg.develop(PackageSpec(path = "/Users/willpope/.julia/dev/BellmanPDEs"))
using BellmanPDEs

include("shield_functions.jl")
include("shield_utils.jl")

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
veh = define_vehicle(wheelbase, body_dims, origin_to_cent, phi_max, v_max)

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


# define human positions and velocity
nearby_human_positions = [[2.0, 7.5],
                        [4.5, 6.0],
                        [4.3, 8.4],
                        [4.6, 3.2],
                        [2.8, 5.1],
                        [4.8, 4.6]]

v_human = 1.0

Dt_plan = 0.5

x_k = SVector(2.0, 2.0, pi*1/3, 1.0)
ia_k = 5
Dt_obs_to_k1 = 0.12
   
# ia_k1_safe_set = shield_action_set(x_k, ia_k, nearby_human_positions, Dt_obs_to_k1, Dt_plan, get_actions, veh)

# @btime shield_action_set($x_k, $ia_k, $nearby_human_positions, $Dt_obs_to_k1, $Dt_plan, $get_actions, $veh)


x_ih_ks1 = nearby_human_positions[2]
x_ih_obs = x_ih_ks1

ig = 1
Dt = 0.5

Dv_max = 0.5
v_k2_max = 1.5
kd_max = ceil(Int, (0.0 - v_k2_max)/(-Dv_max))

# x_ih_ks2 = propagate_human(x_ih_ks1, ig, Dt, v_human, goal_positions)

F_ih_seq, F_ih_body_seq = generate_F_ih_seq(x_ih_obs, Dt_obs_to_k1, Dt_plan, v_human, goal_positions, kd_max)

# F_all_seq = generate_F_all_seq(nearby_human_positions, Dt_obs_to_k1, Dt_plan, v_human, goal_positions, kd_max)