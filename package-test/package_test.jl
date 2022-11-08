# package_test.jl

using Pkg
Pkg.develop(PackageSpec(path = "/Users/willpope/.julia/dev/BellmanPDEs"))
using BellmanPDEs

using BenchmarkTools
using BSON: @save, @load

# define environment (workspace, obstacles, goal)
ws_width = 20.0
ws_length = 20.0
workspace = VPolygon([[0.0, 0.0], [ws_width, 0.0], [ws_width, ws_length], [0.0, ws_length]])

# # standard circular obstacles
# obstacle_list = [VPolyCircle([2.5, 3.0], 0.5), 
#                 VPolyCircle([3.7, 7.5], 0.5),
#                 VPolyCircle([5.9, 4.8], 0.5)]

# goal = VPolyCircle([1/2*ws_width, ws_length-0.75], 0.75)

# local minima
# obstacle_list = [VPolygon([[2.0, 6.0], [6.0, 6.0], [6.0, 7.0], [2.0, 7.0]]),
#                 VPolygon([[1.0, 4.0], [2.0, 4.0], [2.0, 6.0], [1.0, 6.0]]),
#                 VPolygon([[6.0, 4.0], [7.0, 4.0], [7.0, 6.0], [6.0, 6.0]])]

# 20x20 environment
obstacle_list = [VPolyCircle([5.125, 4.875], 1.125), 
                VPolyCircle([6.5, 15.25], 1.5),
                VPolyCircle([16.25, 11.0], 1.125),
                VPolyCircle([10.0, 9.5], 2.25)]

goal = VPolyCircle([2/3*ws_width, ws_length-1.2], 1.2)

# empty
# obstacle_list = []     

env = define_environment(workspace, obstacle_list, goal)

# define vehicle and dynamics
wheelbase = 0.75
body_dims = [1.0, 0.5]
origin_to_cent = [0.375, 0.0]
phi_max = 0.475
v_max = 2.0
veh = define_vehicle(wheelbase, body_dims, origin_to_cent, phi_max, v_max)

Dt = 0.5

# define state grid
state_space = [[0.0, ws_width], [0.0, ws_length], [-pi, pi], [0.0, v_max]]
dx_sizes = [1/2, 1/2, deg2rad(22.5), 1/3]
angle_wrap = [false, false, true, false]
sg = define_state_grid(state_space, dx_sizes, angle_wrap)

# solver params
Dval_tol = 0.01
max_solve_steps = 500

# planner params
max_plan_steps = 2e3

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

    actions = SVector{21, SVector{2, Float64}}(
        (-phi_lim_n, -Dv_lim),        # Dv = -Dv
        (-2/3*phi_lim_n, -Dv_lim),
        (-1/3*phi_lim_n, -Dv_lim),
        (0.0, -Dv_lim),
        (1/3*phi_lim_n, -Dv_lim),
        (2/3*phi_lim_n, -Dv_lim),
        (phi_lim_n, -Dv_lim),

        (-phi_lim, 0.0),       # Dv = 0.0
        (-2/3*phi_lim, 0.0),
        (-1/3*phi_lim, 0.0),
        (0.0, 0.0),
        (1/3*phi_lim, 0.0),
        (2/3*phi_lim, 0.0),
        (phi_lim, 0.0),

        (-phi_lim_p, Dv_lim),        # Dv = +Dv
        (-2/3*phi_lim_p, Dv_lim),
        (-1/3*phi_lim_p, Dv_lim),
        (0.0, Dv_lim),
        (1/3*phi_lim_p, Dv_lim),
        (2/3*phi_lim_p, Dv_lim),
        (phi_lim_p, Dv_lim))

    ia_set = collect(1:length(actions))

    return actions, ia_set
end

# define reactive controller constraints

# define standard HJB reward function
function get_HJB_reward(x, a, Dt, veh)
    reward_x_a = 0
    
    # time penalty
    reward_x_a += -Dt

    return reward_x_a
end

# define POMDP reward function
function get_POMDP_reward(x, a, Dt, veh)
    reward_x_a = 0

    # low speed penalty
    reward_x_a += -(veh.v_max - abs(x[4] + a[2]))/veh.v_max

    # time penalty
    reward_x_a += -1.0

    return reward_x_a
end

# MAIN ---
path = "/Users/willpope/Desktop/Research/marmot-algs/package-test"

solve_flag = 0
plot_value = 0

plan_flag = 1
plot_path = 1

get_reward = get_POMDP_reward

# solver
if solve_flag == 1
    @time q_value_array, value_array = solve_HJB_PDE(get_actions, get_reward, Dt, env, veh, sg, Dval_tol, max_solve_steps)
    @save path * "/bson/HJB_solution.bson" q_value_array value_array
else
    @load path * "/bson/HJB_solution.bson" q_value_array value_array
end

# TO-DO: see how HJB value compares with actual accumulated reward

# planner
if plan_flag == 1
    # x_0_list = [SVector(2.0, 1.5, deg2rad(15), 0.0),
    #             SVector(3.0, 11.5, deg2rad(105), 2.0),
    #             SVector(14.0, 4.5, deg2rad(-60), 0.5),
    #             SVector(10.0, 1.5, deg2rad(90), 1.0)]

    # label_list = ["", "", "", ""]

    # x_path_HJB_1, x_subpath_HJB_1, _, val_path_HJB_1 = plan_path(x_0_list[1], HJB_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    # x_path_HJB_2, x_subpath_HJB_2, _, val_path_HJB_2 = plan_path(x_0_list[2], HJB_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    # x_path_HJB_3, x_subpath_HJB_3, _, val_path_HJB_3 = plan_path(x_0_list[3], HJB_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    # x_path_HJB_4, x_subpath_HJB_4, _, val_path_HJB_4 = plan_path(x_0_list[4], HJB_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)

    # path_list = [x_path_HJB_1, x_path_HJB_2, x_path_HJB_3, x_path_HJB_4]
    # subpath_list = [x_subpath_HJB_1, x_subpath_HJB_2, x_subpath_HJB_3, x_subpath_HJB_4]
    # val_path = [val_path_HJB_1, val_path_HJB_2, val_path_HJB_3, val_path_HJB_4]

    
    # x_0 = SVector(2.0, 1.5, deg2rad(0), 0.0)
    x_0 = SVector(1.8956199596166858, 2.076323753226284, 2.4061404910428497, 1.0)
    label_list = ["Optimal", "Approx Optimal", "Reactive", "Approx Reactive"]
    x_path_HJB, x_subpath_HJB, _, val_path_HJB = plan_path(x_0, HJB_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    x_path_aHJB, x_subpath_aHJB, _, val_path_aHJB = plan_path(x_0, approx_HJB_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    # x_path_RC, x_subpath_RC, _, val_path_RC = plan_path(x_0, reactive_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    # x_path_aRC, x_subpath_aRC, _, val_path_aRC = plan_path(x_0, approx_reactive_policy, get_actions, get_reward, Dt, q_value_array, value_array, env, veh, sg, max_plan_steps)
    path_list = [x_path_HJB, x_path_aHJB]#, x_path_RC, x_path_aRC]
    subpath_list = [x_subpath_HJB, x_subpath_aHJB]#, x_subpath_RC, x_subpath_aRC]
    val_path = [val_path_HJB, val_path_aHJB]#, val_path_RC, val_path_aRC]
end

# plotting
if plot_value == 1
    heatmap_clim = -15
    plot_HJB_value(value_array, env, veh, sg, heatmap_clim)
end

if plot_path == 1
    # plot_path_value([val_path_HJB, val_path_aHJB, val_path_RC, val_path_aRC], Dt)

    linez_clim = 2.5
    plot_HJB_path(path_list, subpath_list, env, veh, linez_clim, label_list)
end

# , x_path_RC, x_path_aRC
# , x_subpath_RC, x_subpath_aRC

# benchmarking
#=
baseline performance:
- discrete_time_EoM:    30.296 ns   (0 allocations)
- propagate_state:      494.078 ns  (8 allocations)
- interp_value:         349.815 ns  (5 allocations)
- optimize_action:      19.243 us   (224 allocations) [21 actions]
- get_actions:          1.044 us    (22 allocations)
- HJB_policy:           20.908 us   (256 allocation)
- reactive_policy:       8.268 us    (110 allocations)

improved performance:
- discrete_time_EoM:    28.974 ns   (0 allocations)
- propagate_state:      163.798 ns  (1 allocation) [think in MVector]
- interp_value:         349.815 ns  (5 allocations) [all from GridInterpolations.jl]
- optimize_action:      9.603 us    (78 allocations) [21 actions] [basically as good as possible]
- get_actions:          49.795 ns   (0 allocations)
- HJB_policy:           10.040 us   (85 allocations) [21 actions]
- approx_HJB_policy:    1.388 us    (23 allocations) [21 actions]
- reactive_policy:      4.285 us    (42 allocations) [21->7 actions]
- approx_policy:        1.430 us    (19 allocations) [21->1 actions, 2^4 neighbors]

structure:
- HJB_policy:                           12.340 us   (110 allocations)
    - get_actions (x1):                 49.795 ns   (0 allocations)
    - optimize_action (x1):             11.844 us   (98 allocations)
        - propagate_state (x21):        163.798 ns  (1 allocation)
            - discrete_time_EoM (x4):   28.974 ns   (0 allocations)
        - interp_value (x21):           349.815 ns  (5 allocations)
=#

# x_k = SVector(8.1, 3.57, deg2rad(35), 1.5)

# actions, ia_set = get_actions(x_k, Dt, veh)

# # a_k = actions[7]
# Dv_RC = -0.5


# debugging issue with LazySets
# veh_body = state_to_body_circle(x_k, veh)

# works fine
# in_target_set(x_k, env, veh)

# works fine
# in_workspace(x_k, env, veh)

# ISSUE: in_obstacle_set() causes Julia to crash
#   - comes from isempty() collision check line
#   - can try isdisjoint() as alternative function
#   - wonder if it's a type issue?

# in_obstacle_set(x_k, env, veh)

# using LazySets
# isempty(intersection(env.obstacle_list[1], veh_body))


# @btime discrete_time_EoM($x_k, $a_k, $Dt, $veh)

# x_k1, x_k1_subpath = propagate_state(x_k, a_k, Dt, veh)
# @btime propagate_state($x_k, $a_k, $Dt, $veh)

# val_x = interpolate(sg.state_grid, value_array, x_k)
# @btime interpolate($sg.state_grid, $value_array, $x_k)

# val_x = interp_value(x_k, value_array, sg)
# @btime interp_value($x_k, $value_array, $sg)

# qvals_x, val_x, ia_opt = optimize_action(x_k, ia_set, actions, get_reward, Dt, value_array, veh, sg)
# @btime optimize_action($x_k, $ia_set, $actions, $get_reward, $Dt, $value_array, $veh, $sg)

# @btime get_actions($x_k, $Dt, $veh)

# @btime HJB_policy($x_k, $Dv_RC, $get_actions, $get_reward, $Dt, $q_value_array, $value_array, $veh, $sg)

# @btime approx_HJB_policy($x_k, $Dv_RC, $get_actions, $get_reward, $Dt, $q_value_array, $value_array, $veh, $sg)

# @btime reactive_policy($x_k, $Dv_RC, $get_actions, $get_reward, $Dt, $value_array, $veh, $sg)   

# approx_reactive_policy(x_k, Dv_RC, get_actions, get_reward, Dt, q_value_array, value_array, veh, sg)
# @btime approx_reactive_policy($x_k, $Dv_RC, $get_actions, $get_reward, $Dt, $q_value_array, $value_array, $veh, $sg)