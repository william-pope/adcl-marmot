# Hamilton-Jacobi-Bellman demonstration
using Plots
using BSON
using BenchmarkTools
using ProfileView

include("HJB_generator_functions.jl")
include("HJB_planner_functions.jl")
include("HJB_utils.jl")
include("dynamics_models.jl")
include("HJB_plotting.jl")

#Returns an environment object
function define_environment()
    # define workspace
    # W = [[0.0 0.0];
    #     [5.518 0.0];
    #     [5.518 11.036];
    #     [0.0 11.036]]

    W = [[0.0 0.0];
        [8.0 0.0];
        [8.0 12.0];
        [0.0 12.0]]

    # define target set
    # T_xy = [[3.15 9.8];
    #         [4.35 9.8];
    #         [4.35 11.0];
    #         [3.15 11.0]]

    T_xy = [[3.25 10.5];
            [4.75 10.5];
            [4.75 12.0];
            [3.25 12.0]]

    T_theta = [[-pi, pi]]

    # define obstacles
    # - circular obstacles defined as [x, y, r], converted to polygon overapproximation
    # OC1 = circle_to_polygon([1.5, 3.0, 0.6])
    # OC2 = circle_to_polygon([4.0, 5.0, 0.6])
    # OC3 = circle_to_polygon([2.8, 8.0, 0.6])

    OC1 = circle_to_polygon([2.5, 3.0, 0.6])
    OC2 = circle_to_polygon([6.0, 5.0, 0.6])
    OC3 = circle_to_polygon([3.8, 8.0, 0.6])

    O_vec = [OC1, OC2, OC3]
    # O_vec = []

    # initialize state grid
    h_xy = 0.125
    h_theta = deg2rad(10)

    env = Environment(h_xy,
                    h_theta,
                    minimum(W[:,1]) : h_xy : maximum(W[:,1]),
                    minimum(W[:,2]) : h_xy : maximum(W[:,2]),
                    -pi : h_theta : pi,
                    W,
                    T_xy,
                    T_theta,
                    O_vec)

    return env
end

#Returns a vehicle object and actions for the vehicle
function define_vehicle()

    # vehicle parameters
    marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)
    unit_car = Vehicle(1.0, 0.5, 0.5, 0.5, 0.75, 0.375, 0.125)
    veh = unit_car

    # Reeds-Shepp car
    # actions = [[a_v,a_phi] for a_v in [-veh.u_vb_max, veh.u_vf_max], a_phi in [-veh.u_phi_max, 0.0, veh.u_phi_max]]

    # Dubins car
    actions = [[a_v, a_phi] for a_v in [veh.u_vf_max], a_phi in [-veh.u_phi_max, -1/2*veh.u_phi_max, 0.0, 1/2*veh.u_phi_max, veh.u_phi_max]]
    actions = reshape(actions, (length(actions),1))
    sort!(actions, dims=1)
    return veh,actions
end

function get_HJB_value_function(env,vehicle,actions,EoM,solve_HJB_flag,plot_growth_flag,max_steps)
    if solve_HJB_flag == true
        du_tol = 0.005
        dt_solve_HJB = 0.5
        value_array, target_array, obstacle_array = solve_HJB_PDE(actions, dt_solve_HJB, du_tol, max_steps, env, vehicle, EoM, plot_growth_flag)
        HJB_dict = Dict(:value_array => value_array, :target_array => target_array, :obstacle_array => obstacle_array,   :env => env, :veh => vehicle)
        bson("./HJB-planner/bson/HJB_dict.bson", HJB_dict)
        N_grid = size(env.x_grid,1) * size(env.y_grid,1) * size(env.theta_grid,1)
        println("total grid nodes = ", N_grid)
        return value_array, target_array, obstacle_array
    else
        HJB_dict = BSON.load("./bson/HJB_dict.bson")
        return HJB_dict[:value_array],HJB_dict[:target_array],HJB_dict[:obstacle_array]
    end
end

function get_HJB_path(starting_point_list, env, vehicle, actions, value_array, obstacle_array, EoM, max_steps)
    dt_generate_HJB_path = 0.5
    if generate_HJB_path_flag == true
        x_path_list = []
        for starting_point in starting_point_list
            println("value at x_0 = ", interp_value(starting_point, value_array, env))
            x_path, u_path, step = plan_HJB_path(starting_point, actions, dt_generate_HJB_path, value_array,
                                                                    obstacle_array, max_steps, EoM, env, vehicle)
            push!(x_path_list, x_path)
            path_time = step*dt_generate_HJB_path
            println("path execution time: ", path_time, " sec")
        end
    end
    return x_path_list
end

#Global Parameters
algs_path_mac = "./"
algs_path_nuc = "/home/adcl/Documents/marmot-algs/"
algs_path = algs_path_mac
solve_HJB_flag = true
generate_HJB_path_flag = true
plot_growth_flag = true
plot_result_flag = true
EoM = dubins_car_EoM
max_steps = 5e3
HJB_path_list = []

# Environment and the vehilcle tosolve the HJB equation for
env = define_environment()
veh,actions = define_vehicle()
value_array,target_array,obstacle_array = get_HJB_value_function(env,veh,actions,EoM,solve_HJB_flag,plot_growth_flag,max_steps)

# Generate paths to goal
starting_point_list = [[6.3, 1.0, deg2rad(90)],
            [4.8, 3.0, deg2rad(105)],
            [3.0, 1.2, deg2rad(80)],
            [1.7, 1.5, deg2rad(165)]]

if generate_HJB_path_flag
    HJB_path_list = get_HJB_path(starting_point_list, env, veh, actions, value_array, obstacle_array, EoM, max_steps )
end

if plot_result_flag
    plot_HJB_result(value_array, HJB_path_list, env, veh)
end

# @btime HJB_action(x_0, value_array, A, obstacle_array, env, veh)

# ProfileView.@profview for _ in 1:1000
#     HJB_action(x_0, value_array, A, obstacle_array, env, veh)
# end



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
