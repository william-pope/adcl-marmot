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
plot_growth_flag = false
plot_value_flag = true
heatmap_clim = 15

dt_solve = 0.25
dval_tol = 0.5
max_solve_steps = 1e3

# planner params
plan_HJB_flag = false

dt_plan = 0.25
max_plan_steps = 2e3



# 2) DEFINITIONS --- --- ---
# define environment (workspace, obstacles, goal)
workspace = VPolygon([[0.0, 0.0], [5.5, 0.0], [5.5, 11.0], [0.0, 11.0]])
obstacle_list = [circle2vpolygon([1.5, 3.0], 0.5), 
                circle2vpolygon([2.7, 7.5], 0.5),
                circle2vpolygon([3.9, 4.8], 0.5)]
goal = VPolygon([[2.25, 10.0], [3.25, 10.0], [3.25, 11.0], [2.25, 11.0]])
env = define_environment(workspace, obstacle_list, goal)

# define vehicle and dynamics
veh_name = "marmot"
veh = define_vehicle(veh_name)

EoM = bicycle_3d_EoM

# define state grid
state_space = [[0.0, 5.5], [0.0, 11.0], [-pi, pi]]
dx_sizes = [0.25, 0.25, deg2rad(10)]
angle_wrap = [false, false, true, false]
sg = define_state_grid(state_space, dx_sizes, angle_wrap)

# define action set
action_space = [[1.0], [-0.475, 0.475]]
du_num_steps = [1, 3]
action_grid = define_action_grid(action_space, du_num_steps)

# define initial states for paths
x_0_list = [[2.3, 1.0, deg2rad(90)],
            [4.8, 3.0, deg2rad(105)],
            [3.0, 1.2, deg2rad(80)],
            [1.7, 1.5, deg2rad(165)]]


# ProfileView.@profview 
# SPEED
#   - 
# @btime initialize_value_array(sg, env, veh)

# TO-DO: 
#   - clean up interpolation
#       - all states in grid neighboring an obstacle node should also be treated as an obstacle node
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

# STATUS (09/01/2022):
#   - solver works with obstacles, solutions look correct
#   - need to work on action space definition and handling

# STATUS (09/22/2022):
#   - solver looks good for 4d_v, some issues with convergence

# TO-DO:
#   ( ) make plotting work

# ISSUE: convergence stalls at dval=0.25 for 4d_v EoM
#   - wonder if related to grid size? or time step?

# 3) MAIN --- --- ---
println("\nstart --- --- ---")

if solve_HJB_flag == true    
    value_array = solve_HJB_PDE(env, veh, EoM, sg, action_grid, dt_solve, dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)

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


# plan paths to goal
if plan_HJB_flag == true
    x_path_list = []

    for x_0 in x_0_list
        # println("value at x_0 = ", interp_value(x_0, value_array, env))

        x_path, u_path, step = plan_HJB_path(x_0, actions, dt_plan, value_array, obstacle_array, max_steps, EoM, env, veh)
        push!(x_path_list, x_path)

        # path_time = step*dt_plan
        # println("path execution time: ", path_time, " sec")
    end
end
    

# # 4) PLOTS --- --- ---

if plot_value_flag == true
    plot_HJB_value(value_array, heatmap_clim, env, veh, sg)
end

if plan_HJB_flag == true
    plot_HJB_path(x_path_list)
end

# # @btime HJB_action(x_0, value_array, A, obstacle_array, env, veh)

# # ProfileView.@profview for _ in 1:1000
# #     HJB_action(x_0, value_array, A, obstacle_array, env, veh)
# # end