# package_test.jl

using Pkg
Pkg.develop(PackageSpec(path = "/Users/willpope/.julia/dev/BellmanPDEs"))
using BellmanPDEs

using BenchmarkTools
using BSON: @save, @load

# define environment (workspace, obstacles, goal)
ws_width = 8.0
ws_length = 12.0
workspace = VPolygon([[0.0, 0.0], [ws_width, 0.0], [ws_width, ws_length], [0.0, ws_length]])
obstacle_list = [VPolyCircle([2.5, 3.0], 0.5), 
                VPolyCircle([3.7, 7.5], 0.5),
                VPolyCircle([5.9, 4.8], 0.5)]
# obstacle_list = []     
goal = VPolyCircle([1/2*ws_width, ws_length-0.75], 0.75)
env = define_environment(workspace, obstacle_list, goal)

# define vehicle and dynamics
wheelbase = 0.324
body_dims = [0.5207, 0.2762]
origin_to_cent = [0.1715, 0.0]
veh = define_vehicle(wheelbase, body_dims, origin_to_cent)

Dt = 0.5

# define state grid
state_space = [[0.0, ws_width], [0.0, ws_length], [-pi, pi], [0.0, 3.0]]
dx_sizes = [0.2, 0.2, deg2rad(15), 1/3]
angle_wrap = [false, false, true, false]
sg = define_state_grid(state_space, dx_sizes, angle_wrap)

# solver params
Dval_tol = 0.01
max_solve_steps = 100

# planner params
max_plan_steps = 2e3

# define action set
function get_actions(x, Dt, veh)
    # set change in velocity (Dv) limit
    Dv_lim = 0.5

    # set steering angle (phi) limit
    phi_max = 0.475
    Dtheta_lim = deg2rad(45)

    v = x[4]
    vp = v + Dv_lim
    vn = v - Dv_lim

    phi_lim = atan(Dtheta_lim * 1/Dt * 1/abs(v) * veh.l)
    phi_lim = clamp(phi_lim, 0.0, phi_max)

    phi_lim_p = atan(Dtheta_lim * 1/Dt * 1/abs(vp) * veh.l)
    phi_lim_p = clamp(phi_lim_p, 0.0, phi_max)

    phi_lim_n = atan(Dtheta_lim * 1/Dt * 1/abs(vn) * veh.l)
    phi_lim_n = clamp(phi_lim_n, 0.0, phi_max)

    actions = SVector{21, SVector{2, Float64}}(
        (-phi_lim, 0.0),       # Dv = 0.0
        (-2/3*phi_lim, 0.0),
        (-1/3*phi_lim, 0.0),
        (0.0, 0.0),
        (1/3*phi_lim, 0.0),
        (2/3*phi_lim, 0.0),
        (phi_lim, 0.0),

        (-phi_lim_n, -Dv_lim),        # Dv = +Dv
        (-2/3*phi_lim_n, -Dv_lim),
        (-1/3*phi_lim_n, -Dv_lim),
        (0.0, -Dv_lim),
        (1/3*phi_lim_n, -Dv_lim),
        (2/3*phi_lim_n, -Dv_lim),
        (phi_lim_n, -Dv_lim),

        (-phi_lim_p, Dv_lim),        # Dv = -Dv
        (-2/3*phi_lim_p, Dv_lim),
        (-1/3*phi_lim_p, Dv_lim),
        (0.0, Dv_lim),
        (1/3*phi_lim_p, Dv_lim),
        (2/3*phi_lim_p, Dv_lim),
        (phi_lim_p, Dv_lim))

    return actions
end

# define cost function
function get_cost(x, a, Dt)
    cost_k = Dt

    return cost_k
end

# MAIN ---
path = "/Users/willpope/Desktop/Research/marmot-algs/package-test"

solve_flag = 1
plot_value = 1

plan_flag = 1
plot_path = 1

# solver
if solve_flag == 1
    @time value_array, a_ind_opt_array = solve_HJB_PDE(get_actions, get_cost, Dt, env, veh, sg, Dval_tol, max_solve_steps)
    @save path * "/bson/HJB_solution.bson" value_array a_ind_opt_array
else
    @load path * "/bson/HJB_solution.bson" value_array a_ind_opt_array
end

# planner
if plan_flag == 1
    x_0 = SVector(2.3, 1.0, deg2rad(15), 2.0)
    x_path, x_subpath, a_path = plan_HJB_path(x_0, get_actions, get_cost, Dt, value_array, a_ind_opt_array, env, veh, sg, max_plan_steps)
end

# plotting
if plot_value == 1
    heatmap_clim = 10
    plot_HJB_value(value_array, heatmap_clim, env, veh, sg)
end

if plot_path == 1
    plot_HJB_path([x_path], [x_subpath], env, veh)
end

# benchmarking
#=
baseline performance:
- discrete_time_EoM:    30.296 ns   (0 allocations)
- propagate_state:      494.078 ns  (8 allocations)
- interp_value:         349.815 ns  (5 allocations)
- optimize_action:      19.243 us   (224 allocations) [21 actions]
- get_actions:          1.044 us    (22 allocations)
- HJB_policy:           20.908 us   (256 allocation)
- rollout_policy:       8.268 us    (110 allocations)

improved performance:
- discrete_time_EoM:    28.974 ns   (0 allocations)
- propagate_state:      163.798 ns  (1 allocation) [think in MVector]
- interp_value:         349.815 ns  (5 allocations) [all from GridInterpolations.jl]
- optimize_action:      11.844 us   (98 allocations) [21 actions] [basically as good as possible]
- get_actions:          49.795 ns   (0 allocations)
- HJB_policy:           12.340 us   (110 allocations) [21 actions]
- rollout_policy:       4.285 us    (42 allocations)

structure:
- HJB_policy:                           12.340 us   (110 allocations)
    - get_actions (x1):                      49.795 ns   (0 allocations)
    - optimize_action (x1):                  11.844 us   (98 allocations)
        - propagate_state (x21):        163.798 ns  (1 allocation)
            - discrete_time_EoM (x4):   28.974 ns   (0 allocations)
        - interp_value (x21):           349.815 ns  (5 allocations)
=#

# actions = get_actions(x_k, Dt, veh)
# a_ind_array = collect(1:length(actions))

# x_k = SVector{4, Float64}(2.3, 1.0, deg2rad(15), 2.0)
# a_k = actions[16]
# Dv_RC = 0.5

# @btime discrete_time_EoM($x_k, $a_k, $Dt, $veh)

# x_k1, x_k1_subpath = propagate_state(x_k, a_k, Dt, veh)
# @btime propagate_state($x_k, $a_k, $Dt, $veh)

# val_x = interpolate(sg.state_grid, value_array, x_k)
# @btime interpolate($sg.state_grid, $value_array, $x_k)

# val_x = interp_value(x_k, value_array, sg)
# @btime interp_value($x_k, $value_array, $sg)

# a_ind_opt, val_x = optimize_action(x_k, a_ind_array, actions, get_cost, Dt, value_array, veh, sg)
# @btime optimize_action($x_k, $a_ind_array, $actions, $get_cost, $Dt, $value_array, $veh, $sg)

# @btime get_actions($x_k, $Dt, $veh)

# @btime HJB_policy($x_k, $get_actions, $get_cost, $Dt, $value_array, $veh, $sg)

# @btime rollout_policy($x_k, $Dv_RC, $get_actions, $get_cost, $Dt, $value_array, $veh, $sg)   