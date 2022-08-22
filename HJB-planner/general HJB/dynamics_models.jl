# dynamics_models.jl

using DomainSets

struct Environment
    workspace::Rectangle    # region
    obstacle_list::Array{Rectangle}     # list of regions
    goal::Rectangle     # region
end

# ISSUE: how will this struct work for non-bicycle model vehicles?
struct Vehicle
    v_range::Array{Float64}     # [m/s]
    phi_range::Array{Float64}   # [rad]
    a_range::Array{Float64}     # [m/s^2]
    xi_range::Array{Float64}    # [rad/s]
    wheelbase::Float64      # wheelbase [m]
    body_length::Float64    # length [m]
    body_width::Float64     # width [m]
    axis_to_back_edge::Float64  # rear bumber to rotation axis [m]
end

struct StateGrid
    grid_array::Array{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
    dx_array::Array{Float64}
end

function gen_environment(env_name)
    if env_name == "aspen_empty"
        workspace = Rectangle((0..5.5), (0..11.0))
        obstacle_list = []
        goal = Rectangle((2.25..3.25), (10.0..11.0))
    end

    env = Environment(workspace, obstacle_list, goal)
    return env
end

function gen_vehicle(veh_name)
    if veh_name == "marmot"
        v_range = [0.0, 1.5]
        phi_range = [-0.475, 0.475]
        a_range = [-1.0, 1.0]
        xi_range = [-3.0, 3.0]
        wheelbase = 0.324
        body_length = 0.5207
        body_width = 0.2762
        axis_to_back_edge = 0.0889

    elseif veh_name == "unit_car"
        v_range = [0.0, 1.0]
        phi_range = [-0.5, 0.5]
        a_range = [-1.0, 1.0]
        xi_range = [-5.0, 5.0]
        wheelbase = 0.5
        body_length = 0.75
        body_width = 0.375
        axis_to_back_edge = 0.125
    end

    veh = Vehicle(v_range, phi_range, a_range, xi_range, wheelbase, body_length, body_width, axis_to_back_edge)
    return veh
end

function gen_state_grid(env, veh, EoM)
    dx_array = [0.25, 0.25, deg2rad(10), 0.1, 0.2]

    x_grid = env.workspace.a[1] : dx[1] : env.workspace.b[1]
    y_grid = env.workspace.a[2] : dx[2] : env.workspace.b[2]
    theta_grid = -pi : dx[3] : pi
    v_grid = veh.v_range[1] : ds[4] : veh.v_range[2]
    phi_grid = veh.phi_range[1] : ds[5] : veh.phi_range[2]

    if EoM == bicycle_3d_EoM
        grid_array = [x_grid, y_grid, theta_grid]
        dx_array = [dx_array[1], dx_array[2], dx_array[3]]

    elseif EoM == bicycle_4d_v_EoM
        grid_array = [x_grid, y_grid, theta_grid, v_grid]
        dx_array = [dx_array[1], dx_array[2], dx_array[3], dx_array[4]]

    elseif EoM == bicycle_4d_phi_EoM
        grid_array = [x_grid, y_grid, theta_grid, phi_grid]
        dx_array = [dx_array[1], dx_array[2], dx_array[3], dx_array[5]]

    elseif EoM == bicycle_5d_EoM
        grid_array = [x_grid, y_grid, theta_grid, v_grid, phi_grid]
        dx_array = [dx_array[1], dx_array[2], dx_array[3], dx_array[4], dx_array[5]]
    end

    state_grid = StateGrid(grid_array, dx_array)
    return state_grid
end

# NOTE: may be able to generalize action/input more, just allowing user to specifiy a discrete (or continuous) set
function gen_actions(veh, EoM)
    v_steps = 2
    phi_steps = 5
    a_steps = 2
    xi_steps = 5

    v_set = collect(range(veh.v_range[1], veh.v_range[2], v_steps))
    phi_set = collect(range(veh.phi_range[1], veh.phi_range[2], phi_steps))
    a_set = collect(range(veh.a_range[1], veh.a_range[2], a_steps))
    xi_set = collect(range(veh.xi_range[1], veh.xi_range[2], xi_steps))

    # TO-DO: need to make these into a list of all possible combinations (need list comprehension syntax)
    if EoM == bicycle_3d_EoM
        actions = [v_set, phi_set]
    elseif EoM == bicycle_4d_v_EoM
        actions = [a_set, phi_set]
    elseif EoM == bicycle_4d_phi_EoM
        actions = [v_set, xi_set]
    elseif EoM == bicycle_5d_EoM
        actions = [a_set, xi_set]
    end

    return actions
end

# equations of motion for 3 DoF kinematic bicycle model
# u = [u_v, u_phi]
function bicycle_3d_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 3)
    
    xdot[1] = u[1]*cos(x[3])    # x
    xdot[2] = u[1]*sin(x[3])    # y
    xdot[3] = u[1]*(1/veh.wheelbase)*tan(u[2])  # theta

    return xdot
end

# equations of motion for 4 DoF velocity kinematic bicycle model
# u = [u_a, u_phi]
function bicycle_4d_v_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 4)
    
    xdot[1] = x[4]*cos(x[3])    # x
    xdot[2] = x[4]*sin(x[3])    # y
    xdot[3] = x[4]*(1/veh.wheelbase)*tan(u[2])  # theta
    xdot[4] = u[1]  # v

    return xdot
end

# equations of motion for 4 DoF steering kinematic bicycle model
# u = [u_v, u_xi]
function bicycle_4d_phi_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 4)
    
    xdot[1] = u[1]*cos(x[3])    # x
    xdot[2] = u[1]*sin(x[3])    # y
    xdot[3] = u[1]*(1/veh.wheelbase)*tan(x[4])  # theta
    xdot[4] = u[2]  # phi

    return xdot
end

# equations of motion for 5 DoF kinematic bicycle model
# u = [u_a, u_xi]
function bicycle_5d_EoM(x::Vector{Float64}, u::Vector{Float64}, veh::Vehicle)
    xdot = zeros(Float64, 5)
    
    xdot[1] = x[4]*cos(x[3])    # x
    xdot[2] = x[4]*sin(x[3])    # y
    xdot[3] = x[4]*(1/veh.wheelbase)*tan(x[5])  # theta
    xdot[4] = u[1]  # v
    xdot[5] = u[2]  # phi

    return xdot
end