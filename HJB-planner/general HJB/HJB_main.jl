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
solve_HJB_flag = false
plot_growth_flag = false
plot_value_flag = false
heatmap_clim = 15

dt_solve = 0.5
dval_tol = 0.01
max_solve_steps = 120

# planner params
plan_HJB_flag = true

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
x_0_list = [[2.8, 0.8, deg2rad(90)],
            [4.1, 1.7, deg2rad(120)],
            [4.8, 2.5, deg2rad(95)],
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


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

if solve_HJB_flag == true    
    value_array, opt_ia_array, set_array = solve_HJB_PDE(env, veh, EoM, sg, ag, dt_solve, dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)

    @save algs_path*"bson/value_array.bson" value_array
    @save algs_path*"bson/opt_ia_array.bson" opt_ia_array
    @save algs_path*"bson/env.bson" env
    @save algs_path*"bson/veh.bson" veh

    # (?): can EoM function be saved as a bson?
    # @save algs_path*"bson/env.bson" EoM
    # @save algs_path*"bson/veh.bson" sg
    # @save algs_path*"bson/env.bson" actions
else
    @load algs_path*"bson/value_array.bson" value_array
    @load algs_path*"bson/opt_ia_array.bson" opt_ia_array
    @load algs_path*"bson/env.bson" env
    @load algs_path*"bson/veh.bson" veh
end

# HJB_policy 
#   -> 230.114 us (full path planner, 3 actions)
#   -> 430.891 us (full path planner, 9 actions)
#   -> 596.315 us (full path planner, 15 actions)
# fast_policy (no safety check)
#   -> 127.109 us (full path planner, 3 actions)
#   -> 119.857 us (full path planner, 9 actions)
#   -> 121.070 us (full path planner, 15 actions)

# x_0 = x_0_list[1]
# @btime x_path, u_path, step = plan_HJB_path(x_0, dt_plan, value_array, opt_ia_array, max_plan_steps, EoM, env, veh, sg, ag)

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