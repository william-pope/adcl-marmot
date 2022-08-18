arr1 = []
@time push!(arr1, 3.0)
@time push!(arr1, 4.0)
@time push!(arr1, 5.0)
@time push!(arr1, 3.0)
@time push!(arr1, 4.0)
@time push!(arr1, 5.0)

arr2 = zeros(Float64, 6)
@time arr2[1] = 3.0
@time arr2[2] = 4.0
@time arr2[3] = 5.0
@time arr2[4] = 3.0
@time arr2[5] = 4.0
@time arr2[6] = 5.0

display(arr1)
display(arr2)