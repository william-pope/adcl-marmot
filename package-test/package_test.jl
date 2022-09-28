# package_test.jl

include("/Users/willpope/.julia/dev/BellmanPDEs/src/BellmanPDEs.jl")

workspace = BellmanPDEs.VPolygon([[0.0, 0.0], [5.5, 0.0], [5.5, 11.0], [0.0, 11.0]])

# looks like package exists, but need less ugly way to use in scripts until listed online