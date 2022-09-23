# HJB_definition_functions.jl

using LazySets
using GridInterpolations

struct Environment
    workspace::VPolygon    # region
    obstacle_list::Array{Any}     # list of regions
    goal::VPolygon     # region
end

# ISSUE: how will this struct work for non-bicycle model vehicles?
#   - in general sense, need a function to map vehicle state to a LazySets polygon for collision/goal checking
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
    state_grid::RectangleGrid
    angle_wrap_array::Array{Bool}
    ind_gs_array::Array
end

function define_environment(workspace, obstacle_list, goal)
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

# TO-DO: catch if state or action dimensions violate the chosen EoM
function define_state_grid(state_space, dx_sizes, angle_wrap)
    state_iters = [minimum(axis):dx_sizes[i]:maximum(axis) for (i, axis) in enumerate(state_space)]
    state_grid = RectangleGrid(state_iters...)

    # Gauss-Seidel
    gs_iters = [[0,1] for axis in state_space]
    gs_prod = Iterators.product(gs_iters...)
    gs_list = Iterators.map(tpl -> convert(SVector{length(gs_iters), Int}, tpl), gs_prod)

    # for sweep in gs_list, need to define ind_list
    ind_gs_array = []
    for (i_gs, gs) in enumerate(gs_list)

        # for axis in sweep = [0,1,1], reverse ind_iters
        ind_iters = Array{StepRange{Int64, Int64}}(undef, size(state_space,1))
        for (i_ax, ax) in enumerate(gs)
            if gs[i_ax] == 0.0
                # forward
                ind_iters[i_ax] = 1:1:size(state_iters[i_ax],1)
            else
                # reverse
                ind_iters[i_ax] = size(state_iters[i_ax],1):-1:1
            end
        end

        ind_prod = Iterators.product(ind_iters...)
        ind_list = Iterators.map(tpl -> convert(SVector{length(ind_iters), Int}, tpl), ind_prod)

        push!(ind_gs_array, ind_list)
    end

    sg = StateGrid(state_grid, angle_wrap, ind_gs_array)
    return sg
end

function define_action_grid(action_space, du_num_steps)
    action_iters = [range(minimum(axis), maximum(axis), du_num_steps[i]) for (i, axis) in enumerate(action_space)]
    action_grid = RectangleGrid(action_iters...)

    return action_grid
end