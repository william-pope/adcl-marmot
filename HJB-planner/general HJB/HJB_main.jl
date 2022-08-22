# Hamilton-Jacobi-Bellman demonstration

using Plots
using BSON: @save, @load
using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_planner_functions.jl")
include("HJB_utils.jl")
include("dynamics_models.jl")
include("HJB_plotting.jl")

# 2) PARAMETERS --- --- ---

 

# # Reeds-Shepp car
# actions = [[a_v,a_phi] for a_v in [-veh.u_vb_max, veh.u_vf_max], a_phi in [-veh.u_phi_max, 0.0, veh.u_phi_max]]

# Dubins car
actions = [[a_v, a_phi] for a_v in [veh.u_vf_max], a_phi in [-veh.u_phi_max, -1/2*veh.u_phi_max, 0.0, 1/2*veh.u_phi_max, veh.u_phi_max]]

actions = reshape(actions, (length(actions),1))
sort!(actions, dims=1)

# define environment (workspace, obstacles, goal_state)
env_name = "aspen_empty"
env = gen_environment(env_name)

# define vehicle and dynamics
veh_name = "marmot"
veh = gen_vehicle(veh_name)

EoM = bicycle_3d_EoM

# define state grid
state_grid = gen_state_grid(env, veh, EoM)

# define action set
actions = gen_actions(veh, EoM)



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

# 3) MAIN --- --- ---
println("\nstart --- --- ---")

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"

algs_path = algs_path_mac

gen_HJB = true
plan_HJB = true

plot_growth = true
plot_result = true

# generate HJB value function
dt_gen = 0.5

if gen_HJB == true
    du_tol = 0.005
    max_steps = 5e3
    
    value_array, target_array, obstacle_array = solve_HJB_PDE(actions, dt_gen, du_tol, max_steps, env, veh, EoM, plot_growth)

    N_grid = size(env.x_grid,1) * size(env.y_grid,1) * size(env.theta_grid,1)
    println("total grid nodes = ", N_grid)

    @save algs_path*"HJB-planner/bson/value_array.bson" value_array
    @save algs_path*"HJB-planner/bson/target_array.bson" target_array
    @save algs_path*"HJB-planner/bson/obstacle_array.bson" obstacle_array
    @save algs_path*"HJB-planner/bson/env.bson" env
    @save algs_path*"HJB-planner/bson/veh.bson" veh
else
    @load algs_path*"HJB-planner/bson/value_array.bson" value_array
    @load algs_path*"HJB-planner/bson/target_array.bson" target_array
    @load algs_path*"HJB-planner/bson/obstacle_array.bson" obstacle_array
    @load algs_path*"HJB-planner/bson/env.bson" env
    @load algs_path*"HJB-planner/bson/veh.bson" veh
end


# plan paths to goal
x_path_list = []

x_0_list = [[6.3, 1.0, deg2rad(90)],
            [4.8, 3.0, deg2rad(105)],
            [3.0, 1.2, deg2rad(80)],
            [1.7, 1.5, deg2rad(165)]]

dt_plan = 0.5

if plan_HJB == true
    max_steps = 5e3

    for x_0 in x_0_list
        println("value at x_0 = ", interp_value(x_0, value_array, env))

        x_path, u_path, step = plan_HJB_path(x_0, actions, dt_plan, value_array, obstacle_array, max_steps, EoM, env, veh)
        push!(x_path_list, x_path)

        path_time = step*dt_plan
        println("path execution time: ", path_time, " sec")
    end
end
    

# 4) PLOTS --- --- ---

plot_HJB_result(value_array, x_path_list, env, veh)

# @btime HJB_action(x_0, value_array, A, obstacle_array, env, veh)

# ProfileView.@profview for _ in 1:1000
#     HJB_action(x_0, value_array, A, obstacle_array, env, veh)
# end