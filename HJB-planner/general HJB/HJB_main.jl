# Hamilton-Jacobi-Bellman test script

using StaticArrays
using BSON: @save, @load
using BenchmarkTools

include("HJB_definition_functions.jl")
include("HJB_solver_functions.jl")
include("HJB_planner_functions.jl")
include("HJB_utils.jl")
include("HJB_plotting.jl")


# 1) GLOBAL PARAMETERS --- --- ---
algs_path_mac = "./"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"
algs_path = algs_path_mac

# solver params
solve_HJB_flag = 0
Dval_tol = 0.01
max_solve_steps = 100

plot_growth_flag = 0
plot_value_flag = 0
heatmap_clim = 10

# planner params
plan_HJB_flag = 1
max_plan_steps = 2e3


# 2) DEFINITIONS --- --- ---
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

# TO-DO: need to make sure that reactive controller velocity still allows vehicle to reach goal

function get_ro_action_set(x)
    # set change in velocity (Dv) limit
    Dv_lim = 0.5

    # set steering angle (phi) limit
    Dt = 0.5
    l = 0.324
    phi_max = 0.475
    Dtheta_lim = deg2rad(45)

    v = x[4]
    vp = v + Dv_lim
    vn = v - Dv_lim

    phi_lim = atan(Dtheta_lim * 1/Dt * 1/abs(v) * l)
    phi_lim = clamp(phi_lim, 0.0, phi_max)

    phi_lim_p = atan(Dtheta_lim * 1/Dt * 1/abs(vp) * l)
    phi_lim_p = clamp(phi_lim_p, 0.0, phi_max)

    phi_lim_n = atan(Dtheta_lim * 1/Dt * 1/abs(vn) * l)
    phi_lim_n = clamp(phi_lim_n, 0.0, phi_max)

    ro_action_set = [[-phi_lim, 0.0],       # Dv = 0.0
                    [-2/3*phi_lim, 0.0],
                    [-1/3*phi_lim, 0.0],
                    [0.0, 0.0],
                    [1/3*phi_lim, 0.0],
                    [2/3*phi_lim, 0.0],
                    [phi_lim, 0.0],

                    [-phi_lim_n, -Dv_lim],        # Dv = +Dv
                    [-2/3*phi_lim_n, -Dv_lim],
                    [-1/3*phi_lim_n, -Dv_lim],
                    [0.0, -Dv_lim],
                    [1/3*phi_lim_n, -Dv_lim],
                    [2/3*phi_lim_n, -Dv_lim],
                    [phi_lim_n, -Dv_lim],

                    [-phi_lim_p, Dv_lim],        # Dv = -Dv
                    [-2/3*phi_lim_p, Dv_lim],
                    [-1/3*phi_lim_p, Dv_lim],
                    [0.0, Dv_lim],
                    [1/3*phi_lim_p, Dv_lim],
                    [2/3*phi_lim_p, Dv_lim],
                    [phi_lim_p, Dv_lim]]

    return ro_action_set
end

# define cost function
function get_cost(x, a, Dt)
    cost_k = Dt

    return cost_k
end

# define initial states for path planner
# x_0_list = [[2.8, 0.8, deg2rad(90)],
#             [4.1, 1.7, deg2rad(120)],
#             [4.8, 2.5, deg2rad(95)],
#             [1.7, 1.5, deg2rad(165)]]

x_0_list = [SA[5.8, 2.0, deg2rad(-120), 1.0]]


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

if solve_HJB_flag == 1    
    # ProfileView.@profview solve_HJB_PDE(env, veh, sg, Dt, Dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)
    @time value_array, a_ind_opt_array, set_array = solve_HJB_PDE(env, veh, sg, Dt, Dval_tol, max_solve_steps, plot_growth_flag, heatmap_clim)

    @save algs_path*"bson/value_array.bson" value_array
    @save algs_path*"bson/a_ind_opt_array.bson" a_ind_opt_array
else
    @load algs_path*"bson/value_array.bson" value_array
    @load algs_path*"bson/a_ind_opt_array.bson" a_ind_opt_array
end

# plan paths to goal
if plan_HJB_flag == 1
    x_path_list = []
    x_subpath_list = []

    for x_0 in x_0_list
        # println("value at x_0 = ", interp_value(x_0, value_array, env))

        x_path, x_subpath, a_path, step = plan_rollout_path(x_0, Dt, value_array, a_ind_opt_array, max_plan_steps, env, veh, sg)
        x_path_HJB, x_subpath_HJB, _, _ = plan_HJB_path(x_0, Dt, value_array, a_ind_opt_array, max_plan_steps, env, veh, sg)
        
        # display(x_path)
        # display(a_path)

        # for kk in eachindex(x_subpath)
        #     println("kk = ",  kk,  ", x_kk = ", x_subpath[kk])
        # end
        
        push!(x_path_list, x_path)
        push!(x_subpath_list, x_subpath)

        push!(x_path_list, x_path_HJB)
        push!(x_subpath_list, x_subpath_HJB)

        # path_time = step*Dt
        # println("path execution time: ", path_time, " sec")
    end
end
    

# 4) PLOTS --- --- ---
if plot_value_flag == 1
    plot_HJB_value(value_array, heatmap_clim, env, veh, sg)
end

if plan_HJB_flag == 1
    plot_HJB_path(x_path_list, x_subpath_list)
end


# x_k = x_0_list[1]
# Dv_RC = 0.5

# HJB_rollout_policy(x_k, Dv_RC, Dt, value_array, veh, sg)

# (?): what happens when v=0 and RC requests -0.5?
#   - recognizes propagated state v_p=-0.5 as being invalid due to 1e5 value