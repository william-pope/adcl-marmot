# minkowski_test.jl

using LazySets
using Plots

using Pkg
Pkg.develop(PackageSpec(path = "/Users/willpope/.julia/dev/BellmanPDEs"))
using BellmanPDEs

#= NOTES
- MinkowskiSum() (lazy operation) works just fine in here
- minkowski_sum() (concrete operation) also works
- circle needs to be centered at [0.0, 0.0], otherwise msum is shifted
=#

# define polygon
x_k_points = [[1.0, 1.0],
            [2.0, 1.0],
            [2.5, 2.0],
            [1.0, 2.0]]

poly1 = VPolygon(x_k_points)

# define circular polygon
cir1 = VPolyCircle([0.0, 0.0], 0.25)

# take Minkowski sum
# msum1 = MinkowskiSum(poly1, cir1)
msum1 = minkowski_sum(poly1, cir1)


# plotting
p1 = plot(aspect_ratio=:equal)

plot!(p1, poly1)
plot!(p1, cir1)
plot!(p1, msum1)