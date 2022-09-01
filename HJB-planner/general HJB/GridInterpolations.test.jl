# GridInterp test

using GridInterpolations

# rectangular grid can be defined in n-dimensions
# NOTE: grid definition not able to handle backwards arrays
state_grid = RectangleGrid(0:1, 0:1, 0:1) 
value_array = [8., 1., 6., 3., 5., 7., 4., 9.]

x = [1.5, 1, 1]
@show interpolate(state_grid, value_array, x)

# grid data listed as [1,1], [2,1], [3,1], [1,2], [2,2], [3,2], ...
#   - seems difficult to work with
#   - need to access 1-d array based on n-d grid indices
#   - (!): same order as in sg.ind_gs_array[1]
#       - not that useful, because no easy way to look up index within ind_list
#       - should be able to write explicit equation that calculates 1-d index from n-d index



function multi2single_ind(ind_m, sg)
    ind_s = 1
    for d in eachindex(ind_m)
        ind_s += (ind_m[d]-1)*prod(sg.state_grid.cut_counts[1:(d-1)])
    end

    return ind_s
end

ind_m = [3, 2]
lens = [3, 4]

# @show ind_s = multi2single_2(ind_m, lens)