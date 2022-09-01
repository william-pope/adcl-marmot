# Hamilton-Jacobi-Bellman demonstration

using Plots
using DomainSets
using BSON: @save, @load
using BenchmarkTools
using ProfileView

include("HJB_definition_functions.jl")
include("HJB_generator_functions.jl")
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
plot_growth_flag = true
plot_result_flag = true

dt_solve = 0.25
dval_tol = 0.005
max_solve_steps = 1

# planner params
run_HJB_planner_flag = false
HJB_path_list = []

dt_plan = 0.25
max_plan_steps = 2e3


# 2) DEFINITIONS --- --- ---
# define environment (workspace, obstacles, goal)
env_name = "aspen_empty"
env = define_environment(env_name)

# define vehicle and dynamics
veh_name = "marmot"
veh = define_vehicle(veh_name)

EoM = bicycle_3d_EoM

# define state grid
state_space = [[0.0, 5.5], [0.0, 11.0], [-pi/2, pi/2]]
axis_step_sizes = [0.25, 0.25, deg2rad(10)]
angle_wrap = [false, false, true]
sg = define_state_grid(state_space, axis_step_sizes, angle_wrap)

# define action set
actions = define_actions(veh, EoM)
# ISSUE: need to differentiate between Dubins car and Reeds-Shepp car somewhere


# TO-DO: 
#   - clean up interpolation
#       - does x_p need to be collision checked when building value array?
#           - does this add significant runtime to HJB generation?
#           - need to store info somewhere to not repeat checking
#           - determine if this is actually helpful before putting time into it
#       - script works, but not sure if more robust approach needed for obstacle value backups
#   - clean up static types
#   - test methods for faster action search
#   - test differential drive dynamics model
#   - build package out of code
#   - replace matrices with arrays (can this be done?)
#   - make Vehicle struct, handlings more general
#   - add action set definition to Vehicle struct (should this be done?)
#       - would be better if it was one big vehicle package (?) for:
#           - EoM
#           - action set
#           - physical parameters

# STATUS:
#   - initialize_value_array() seems to be working, but haven't checked rigorously
# #   - 5d model takes forever to solve

# STATUS (08/31/22):
#   - initialize() works with new GridInterpolations method (checked free space, edge, goal values)

# 3) MAIN --- --- ---
println("\nstart --- --- ---")

if solve_HJB_flag == true    
    value_array = solve_HJB_PDE(env, veh, EoM, sg, actions, dt_solve, dval_tol, max_solve_steps, plot_growth_flag)

    @save algs_path*"bson/value_array.bson" value_array
    @save algs_path*"bson/env.bson" env
    @save algs_path*"bson/veh.bson" veh

    # (?): can EoM function be saved as a bson?
    # @save algs_path*"bson/env.bson" EoM
    # @save algs_path*"bson/veh.bson" sg
    # @save algs_path*"bson/env.bson" actions
else
    @load algs_path*"bson/value_array.bson" value_array
    @load algs_path*"bson/env.bson" env
    @load algs_path*"bson/veh.bson" veh
end


# # plan paths to goal
x_path_list = []

# x_0_list = [[6.3, 1.0, deg2rad(90)],
#             [4.8, 3.0, deg2rad(105)],
#             [3.0, 1.2, deg2rad(80)],
#             [1.7, 1.5, deg2rad(165)]]

# dt_plan = 0.5

# if plan_HJB == true
#     max_steps = 5e3

#     for x_0 in x_0_list
#         println("value at x_0 = ", interp_value(x_0, value_array, env))

#         x_path, u_path, step = plan_HJB_path(x_0, actions, dt_plan, value_array, obstacle_array, max_steps, EoM, env, veh)
#         push!(x_path_list, x_path)

#         path_time = step*dt_plan
#         println("path execution time: ", path_time, " sec")
#     end
# end
    

# # 4) PLOTS --- --- ---

# plot_HJB_result(value_array, x_path_list, env, veh)

# # @btime HJB_action(x_0, value_array, A, obstacle_array, env, veh)

# # ProfileView.@profview for _ in 1:1000
# #     HJB_action(x_0, value_array, A, obstacle_array, env, veh)
# # end