# package_test.jl

using Pkg
Pkg.develop(PackageSpec(path = "/Users/willpope/.julia/dev/BellmanPDEs"))

using BellmanPDEs
using BenchmarkTools

# define environment (workspace, obstacles, goal)
workspace = VPolygon([[0.0, 0.0], [5.5, 0.0], [5.5, 11.0], [0.0, 11.0]])
obstacle_list = [VPolyCircle([1.5, 3.0], 0.5), 
                VPolyCircle([2.7, 7.5], 0.5),
                VPolyCircle([3.9, 4.8], 0.5)]
goal = VPolygon([[2.125, 9.75], [3.375, 9.75], [3.375, 11.0], [2.125, 11.0]])
env = define_environment(workspace, obstacle_list, goal)

# define vehicle and dynamics
wheelbase = 0.324
body_dims = [0.5207, 0.2762]
origin_to_cent = [0.1715, 0.0]
veh = define_vehicle(wheelbase, body_dims, origin_to_cent)

function bicycle_3d_EoM(x::SVector{3, Float64}, u::SVector{2, Float64}, veh::Vehicle)
    x_dot1 = u[1]*cos(x[3])
    x_dot2 = u[1]*sin(x[3])
    x_dot3 = u[1]*(1/veh.wheelbase)*tan(u[2])

    x_dot = SVector(x_dot1, x_dot2, x_dot3)

    return x_dot
end

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

# set solver params
dt_solve = 0.5
dval_tol = 0.1
max_solve_steps = 100

# solver
value_array, opt_ia_array, set_array = solve_HJB_PDE(env, veh, EoM, sg, ag, dt_solve, dval_tol, max_solve_steps)

# set policy params
dt_plan = 0.5
x = SVector(3.21, 1.5, deg2rad(60))

# policy
a = fast_policy(x, dt_plan, value_array, opt_ia_array, EoM, veh, sg, ag)