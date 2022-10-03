# Hamilton-Jacobi-Bellman demonstration

using StaticArrays
using BSON: @save, @load
using BenchmarkTools

include("HJB_definition_functions.jl")
include("HJB_solver_functions.jl")
include("HJB_planner_functions.jl")
include("HJB_utils.jl")
include("dynamics_models.jl")
include("HJB_plotting.jl")


# 1) GLOBAL PARAMETERS --- --- ---
algs_path_mac = "./"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"
algs_path = algs_path_mac

# solver params
solve_HJB_flag = true
plot_growth_flag = false
plot_value_flag = false
heatmap_clim = 15

dt_solve = 0.5
dval_tol = 0.1
max_solve_steps = 100

# planner params
plan_HJB_flag = false

dt_plan = 0.5
max_plan_steps = 2e3


# 2) DEFINITIONS --- --- ---
# define environment (workspace, obstacles, goal)
workspace = VPolygon([[0.0, 0.0], [5.5, 0.0], [5.5, 11.0], [0.0, 11.0]])
obstacle_list = [circle2vpolygon([1.5, 3.0], 0.5), 
                circle2vpolygon([2.7, 7.5], 0.5),
                circle2vpolygon([3.9, 4.8], 0.5)]
goal = VPolygon([[2.125, 9.75], [3.375, 9.75], [3.375, 11.0], [2.125, 11.0]])
env = define_environment(workspace, obstacle_list, goal)

# define vehicle and dynamics
wheelbase = 0.324
body_dims = [0.5207, 0.2762]
origin_to_cent = [0.1715, 0.0]
veh = define_vehicle(wheelbase, body_dims, origin_to_cent)

EoM = bicycle_3d_EoM

# define state grid
state_space = [[0.0, 5.5], [0.0, 11.0], [-pi, pi]]
dx_sizes = [0.2, 0.2, deg2rad(15)]
angle_wrap = [false, false, true]
sg = define_state_grid(state_space, dx_sizes, angle_wrap)

# define action set
action_space = [[1.0], [-0.475, 0.475]]
du_num_steps = [1, 5]
ag = define_action_grid(action_space, du_num_steps)

# define initial states for paths
# x_0_list = [[2.8, 0.8, deg2rad(90)],
#             [4.1, 1.7, deg2rad(120)],
#             [4.8, 2.5, deg2rad(95)],
#             [1.7, 1.5, deg2rad(165)]]

x_0_list = [SA[4.1, 1.7, deg2rad(120)]]


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

if solve_HJB_flag == true    
    # ProfileView.@profview solve_HJB_PDE(env, veh, EoM, sg, ag, dt_solve, dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)
    @time value_array, opt_ia_array, set_array = solve_HJB_PDE(env, veh, EoM, sg, ag, dt_solve, dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)

    @save algs_path*"bson/value_array.bson" value_array
    @save algs_path*"bson/opt_ia_array.bson" opt_ia_array
else
    @load algs_path*"bson/value_array.bson" value_array
    @load algs_path*"bson/opt_ia_array.bson" opt_ia_array
end

# plan paths to goal
if plan_HJB_flag == true
    x_path_list = []

    for x_0 in x_0_list
        # println("value at x_0 = ", interp_value(x_0, value_array, env))

        x_path, u_path, step = plan_HJB_path(x_0, dt_plan, value_array, opt_ia_array, max_plan_steps, EoM, env, veh, sg, ag)
        display(u_path)
        
        push!(x_path_list, x_path)

        # path_time = step*dt_plan
        # println("path execution time: ", path_time, " sec")
    end
end
    

# 4) PLOTS --- --- ---
if plot_value_flag == true
    plot_HJB_value(value_array, heatmap_clim, env, veh, sg)
end

if plan_HJB_flag == true
    plot_HJB_path(x_path_list)
end

# x_k = SA[4.1, 1.7, deg2rad(120)]

# @btime HJB_policy(x_k, dt_plan, value_array, EoM, veh, sg, ag)
# @btime fast_policy(x_k, dt_plan, value_array, opt_ia_array, EoM, veh, sg, ag)

# ProfileView.@profview for _ in 1:1000
#     fast_policy(x_k, dt_plan, value_array, opt_ia_array, EoM, veh, sg, ag)
# end

# X) PERFORMANCE --- --- ---


# x_dot = [0.0, 0.0, 0.0]
# x = [3.1, 1.2, 0.7*pi]
# a = [1.0, -0.475]

# x_dot = MVector(0.0, 0.0, 0.0)
x = SVector(3.1, 1.2, 0.7*pi)
a = SVector(1.0, -0.475)

# ProfileView.@profview 

# @btime bicycle_3d_EoM(x, a, veh)
# @btime bicycle_3d_EoM!(x_dot, x, a, veh)

# println("\nEoM")
# @btime bicycle_3d_EoM($x, $a, $veh)

# println("\nRK4")
# @btime runge_kutta_4($x, $a, $dt_solve, $EoM, $veh, $sg)

# println("\ninterp")
# @btime interpolate($sg.state_grid, $value_array, $x)
# @btime interp_value($x, $value_array, $sg)

# println("\npolicy")
# @btime HJB_policy($x, $dt_plan, $value_array, $EoM, $veh, $sg, $ag)
# @btime fast_policy($x, $dt_plan, $value_array, $opt_ia_array, $EoM, $veh, $sg, $ag)

# @profview fast_policy(x, dt_plan, value_array, opt_ia_array, EoM, veh, sg, ag)