# HJB_definition_functions.jl

using DomainSets
using LazySets

struct Environment
    workspace::VPolygon    # region
    obstacle_list::Array{VPolygon}     # list of regions
    goal::VPolygon     # region
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
    axis_to_cent_x::Float64
    axis_to_cent_y::Float64
    origin_body::VPolygon
end

struct StateGrid
    grid_array::Array{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
    angle_wrap_array::Array{Bool}
    # step_array::Array{Float64}
    grid_size_array::Array{Int64}
    state_dim::Int64
    # iter_array::Vector{UnitRange{Int64}}
    grid_idx_matrix::Matrix{Int64}
    gs_sweep_matrix::Matrix{Int64}
end

function define_environment(env_name)
    if env_name == "aspen_empty"
        workspace = VPolygon([[0.0, 0.0], [5.5, 0.0], [5.5, 11.0], [0.0, 11.0]])
        obstacle_list = []
        goal = VPolygon([[2.25, 10.0], [3.25, 10.0], [3.25, 11.0], [2.25, 11.0]])
    end

    env = Environment(workspace, obstacle_list, goal)
    return env
end

function define_vehicle(veh_name)
    if veh_name == "marmot"
        v_range = [-0.75, 1.5]
        phi_range = [-0.475, 0.475]
        a_range = [-1.0, 1.0]
        xi_range = [-3.0, 3.0]
        wheelbase = 0.324
        body_length = 0.5207
        body_width = 0.2762
        axis_to_cent_x = 0.1715
        axis_to_cent_y = 0.0

    elseif veh_name == "unit_car"
        v_range = [-0.5, 1.0]
        phi_range = [-0.5, 0.5]
        a_range = [-1.0, 1.0]
        xi_range = [-5.0, 5.0]
        wheelbase = 0.5
        body_length = 0.75
        body_width = 0.375
        axis_to_cent_x = 0.25
        axis_to_cent_y = 0.0
    end

    x0_min = axis_to_cent_x - 1/2*body_length
    x0_max = axis_to_cent_x + 1/2*body_length
    y0_min = axis_to_cent_y - 1/2*body_width
    y0_max = axis_to_cent_y + 1/2*body_width
    origin_body = VPolygon([[x0_min, y0_min], [x0_max, y0_min], [x0_max, y0_max], [x0_min, y0_max]])

    veh = Vehicle(v_range, phi_range, a_range, xi_range, wheelbase, body_length, body_width, axis_to_cent_x, axis_to_cent_y, origin_body)
    return veh
end

# TO-DO: add array that encodes whether grid wraps around unit circle (theta, phi)
# a[:,2] = reverse(a[:,2]) -> use for Gauss-Seidel sweeps
function define_state_grid(env, veh, EoM)
    step_array = [0.25, 0.25, deg2rad(11.25), 0.25, 0.1]

    # define state grid
    x_grid = minimum(getindex.(env.workspace.vertices, 1)) : step_array[1] : maximum(getindex.(env.workspace.vertices, 1))
    y_grid = minimum(getindex.(env.workspace.vertices, 2)) : step_array[2] : maximum(getindex.(env.workspace.vertices, 2))
    theta_grid = -pi : step_array[3] : pi
    v_grid = veh.v_range[1] : step_array[4] : veh.v_range[2]
    phi_grid = veh.phi_range[1] : step_array[5] : veh.phi_range[2]

    grid_array = [x_grid, y_grid, theta_grid, v_grid, phi_grid]
    angle_wrap_array = [false, false, true, false, true]

    grid_size_array = [size(axis,1) for axis in grid_array]
    iter_array = [1:ax_size for ax_size in grid_size_array]

    # modify arrays for chosen EoM
    if EoM == bicycle_3d_EoM
        delete_idx = [4,5]
    elseif EoM == bicycle_4d_v_EoM
        delete_idx = [5]
    elseif EoM == bicycle_4d_phi_EoM
        delete_idx = [4]
    elseif EoM == bicycle_5d_EoM
        delete_idx = []
    end

    deleteat!(grid_array, delete_idx)
    deleteat!(angle_wrap_array, delete_idx)
    deleteat!(step_array, delete_idx)
    deleteat!(grid_size_array, delete_idx)
    deleteat!(iter_array, delete_idx)

    # define grid_idx_matrix
    state_dim = size(grid_array, 1)
    num_nodes = prod(grid_size_array)
    
    arr1 = collect(Iterators.product(iter_array...))
    arr2 = reshape(arr1, (length(arr1),1))
    
    grid_idx_matrix = zeros(Int64, (num_nodes, state_dim))
    
    for row in eachindex(arr2)
        for d in 1:state_dim
            grid_idx_matrix[row, d] = arr2[row][d]
        end
    end

    # define Gauss-Seidel matrix
    gs_array = fill(0:1, state_dim)
    arr1 = collect(Iterators.product(gs_array...))
    arr2 = reshape(arr1, (length(arr1),1))

    gs_sweep_matrix = zeros(Int64, (2^state_dim, state_dim))

    for row in eachindex(arr2)
        for b in 1:state_dim
            gs_sweep_matrix[row, b] = arr2[row][b]
        end
    end

    sg = StateGrid(grid_array, angle_wrap_array, grid_size_array, state_dim, grid_idx_matrix, gs_sweep_matrix)
    return sg
end

# NOTE: should be able to generalize action/input more, just allowing user to specifiy a discrete (or continuous) set
#   - would this be better? probably
function define_actions(veh, EoM)
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
        action_list = [v_set, phi_set]
    elseif EoM == bicycle_4d_v_EoM
        action_list = [a_set, phi_set]
    elseif EoM == bicycle_4d_phi_EoM
        action_list = [v_set, xi_set]
    elseif EoM == bicycle_5d_EoM
        action_list = [a_set, xi_set]
    end

    actions = [[a1, a2] for a1 in action_list[1], a2 in action_list[2]]     # NOTE: probably a better way to do this
    actions = reshape(actions, (length(actions),1))
    sort!(actions, dims=1)

    return actions
end