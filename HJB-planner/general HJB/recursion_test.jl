# recursion test

# x_grid = 0:0.25:0.25
# y_grid = 0:0.5:1
# th_grid = 0:0.5:1

# grid_array = [x_grid, y_grid, th_grid]
# len_array = [size(x_grid,1), size(y_grid,1), size(th_grid,1)]
# idx_array = (1:len_array[1], 1:len_array[2], 1:len_array[3])

# state_dim = size(grid_array,1)
# num_nodes = prod(len_array)

# # direct method
# g = collect(Iterators.product(idx_array...))
# k = reshape(g, (length(g),1))

# m = zeros(Int64, (num_nodes, state_dim))

# for i in eachindex(k)
#     for d in eachindex(grid_array)
#         m[i,d] = k[i][d]
#     end
# end

# display(m)

state_dim = 5

gs_array = fill(0:1, state_dim)
arr1 = collect(Iterators.product(gs_array...))
arr2 = reshape(arr1, (length(arr1),1))

gs_sweep_matrix = zeros(Int64, (2^state_dim, state_dim))

for row in eachindex(arr2)
    for b in 1:state_dim
        gs_sweep_matrix[row, b] = arr2[row][b]
    end
end

display(gs_sweep_matrix)