# state grid test

using StaticArrays
using DomainSets
using GridInterpolations

state_space = [[0.0, 1.0], [0.0, 1.5], [0.0, 2.0]]
state_steps = [0.5, 0.5, 0.5]
angle_wrap = [false, false, true]

state_iters = [minimum(axis):state_steps[i]:maximum(axis) for (i, axis) in enumerate(state_space)]
state_grid = RectangleGrid(state_iters...)

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

println("")
for ind in ind_gs_array[1]
    println(ind, " -> ", state_grid[ind...])
end

# for ind in ind_gs_array[4]
#     println(ind)
# end