# HJB_definition_functions.jl

using DomainSets
using LazySets
using GridInterpolations

struct Environment
    workspace::VPolygon    # region
    obstacle_list::Array{Any}     # list of regions
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

# TO-DO: catch if state or aciton dimensions violate the chosen EoM
# a[:,2] = reverse(a[:,2]) -> use for Gauss-Seidel sweeps
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

# TO-DO: define action space using Range(a,b,n) type
function define_action_grid(action_space, du_nums)
    action_iters = [minimum(axis):du_sizes[i]:maximum(axis) for (i, axis) in enumerate(action_space)]
    action_grid = RectangleGrid(action_iters...)

    return action_grid
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