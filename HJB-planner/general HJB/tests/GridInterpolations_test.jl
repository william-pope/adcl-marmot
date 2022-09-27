# GridInterp test

using GridInterpolations
using StaticArrays

# rectangular grid can be defined in n-dimensions
# NOTE: grid definition not able to handle backwards arrays
state_grid = RectangleGrid(0:1, 0:1, 0:1) 
value_array = [2, 5, 2, 3, 4, 1, 3, 6]

full_set = collect(1:10)

x = [0.18, 0.12, 0.15]
@show interpolate(state_grid, value_array, x)

index, weight = interpolants(state_grid, x)
@show index
@show weight
@show nbr_index_sorted = index[sortperm(weight, rev=true)]     # returns grid indices sorted by distance to state
@show nbr_value_srt_unq = unique(value_array[nbr_index_sorted])

println(" ")
@show leftover_values = setdiff(full_set, value_array)
@show leftover_values_shuf = shuffle(leftover_values)

println(" ")
@show final_order = vcat(nbr_value_srt_unq, leftover_values_shuf)

# @btime nbr_value_sorted = value_array[index[sortperm(weight, rev=true)]]

# permutation vector: elements in vector show what original element is called to that position in permutation vector

# grid data listed as [1,1], [2,1], [3,1], [1,2], [2,2], [3,2], ...
#   - seems difficult to work with
#   - need to access 1-d array based on n-d grid indices
#   - (!): same order as in sg.ind_gs_array[1]
#       - not that useful, because no easy way to look up index within ind_list
#       - should be able to write explicit equation that calculates 1-d index from n-d index

# function multi2single_ind(ind_m, sg)
#     ind_s = 1
#     for d in eachindex(ind_m)
#         ind_s += (ind_m[d]-1)*prod(sg.state_grid.cut_counts[1:(d-1)])
#     end

#     return ind_s
# end

# ind_m = [3, 2]
# lens = [3, 4]

# @show ind_s = multi2single_2(ind_m, lens)

# testing 1-d data array -> n-d data array
# value_array_m = reshape(value_array, state_grid.cut_counts...)