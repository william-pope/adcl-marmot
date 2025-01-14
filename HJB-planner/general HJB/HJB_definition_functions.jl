# HJB_definition_functions.jl

using LazySets
using GridInterpolations
using StaticArrays

struct Environment
    workspace::VPolygon
    obstacle_list::Array{VPolygon}
    goal::VPolygon
end

struct VehicleBody
    wheelbase::Float64              # wheelbase [m]
    body_dims::Array{Float64}       # [length, width] [m]
    origin_to_cent::Array{Float64}  # [x, y] [m]
    origin_body::VPolygon
end

struct StateGrid
    state_grid::RectangleGrid
    state_list_static::Array{Any}
    angle_wrap_array::Array{Bool}
    ind_gs_array::Array
end

# defines environment geoemtry
function define_environment(workspace, obstacle_list, goal)
    env = Environment(workspace, obstacle_list, goal)
    return env
end

# defines vehicle geometry
function define_vehicle(wheelbase, body_dims, origin_to_cent)
    x0_min = origin_to_cent[1] - 1/2*body_dims[1]
    x0_max = origin_to_cent[1] + 1/2*body_dims[1]
    y0_min = origin_to_cent[2] - 1/2*body_dims[2]
    y0_max = origin_to_cent[2] + 1/2*body_dims[2]
    origin_body = VPolygon([[x0_min, y0_min], [x0_max, y0_min], [x0_max, y0_max], [x0_min, y0_max]])

    veh = VehicleBody(wheelbase, body_dims, origin_to_cent, origin_body)
    return veh
end

# discretizes state space
function define_state_grid(state_space, dx_sizes, angle_wrap)
    state_iters = [minimum(axis):dx_sizes[i]:maximum(axis) for (i, axis) in enumerate(state_space)]
    state_grid = RectangleGrid(state_iters...)

    state_list_static = []
    for state in state_grid
        push!(state_list_static, SA[state...])
    end

    # Gauss-Seidel sweeping scheme
    gs_iters = [[0,1] for axis in state_space]
    gs_prod = Iterators.product(gs_iters...)
    gs_list = Iterators.map(tpl -> convert(SVector{length(gs_iters), Int}, tpl), gs_prod)

    # for sweep in gs_list, need to define ind_list
    ind_gs_array = []
    for (i_gs, gs) in enumerate(gs_list)
        ind_iters = Array{StepRange{Int64, Int64}}(undef, size(state_space, 1))
        for (i_ax, ax) in enumerate(gs)
            num_iters = size(state_iters[i_ax], 1)
            if gs[i_ax] == 0.0
                # forward
                ind_iters[i_ax] = 1:1:num_iters
            else
                # reverse
                ind_iters[i_ax] = num_iters:-1:1
            end
        end

        ind_prod = Iterators.product(ind_iters...)
        ind_list = Iterators.map(tpl -> convert(SVector{length(ind_iters), Int}, tpl), ind_prod)

        push!(ind_gs_array, ind_list)
    end

    sg = StateGrid(state_grid, state_list_static, angle_wrap, ind_gs_array)
    return sg
end