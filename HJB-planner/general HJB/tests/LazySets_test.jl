# LazySets.jl test

using LazySets
using Plots
using BenchmarkTools

# NOTE: dimensionality
#   - collision checking (polygon intersection) will be done in 2-d workspace (using [x, y, theta] for vehicle)
#   - in() set checking will be done in full n-dim space to check if state is within n-dim goal region, n-dim state space

function circle2vpolygon(cent_cir, r_cir)
    # number of points used to discretize edge of circle
    pts = 16

    # circle radius is used as midpoint radius for polygon faces (over-approximation)
    r_poly = r_cir/cos(pi/pts)

    theta_rng = range(0, 2*pi, length=pts+1)

    cir_vertices = [[cent_cir[1] + r_poly*cos(theta), cent_cir[2] + r_poly*sin(theta)] for theta in theta_rng]
    
    poly_cir = VPolygon(cir_vertices)
    return poly_cir
end

# define polygons
poly1 = VPolygon([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
poly2 = VPolygon([[0.5, 0.5], [0.9, 0.5], [0.9, 0.9], [0.5, 0.9]])

# define circles
cir1 = circle2vpolygon([1.5, 1.3], 0.5)
cir0 = Ball2([1.5, 1.3], 0.5)

# find intersection of two polygons
int1 = intersection(poly1, poly2)

# test if subset
ss = issubset(poly2, poly1)

# test if disjoint
dj = isdisjoint(poly1, cir1)

# plotting
p1 = plot(poly1, aspect_ratio=:equal)
# plot!(p1, poly2)
plot!(p1, cir1)
plot!(p1, cir0)

display(p1)
display(dj)

# vehicle body transformation function
function state_to_body(x, body)
    # rotate body about origin by theta
    rot_matrix = [cos(x[3]) -sin(x[3]); sin(x[3]) cos(x[3])]
    body = linear_map(rot_matrix, body)

    # translate body from origin by [x, y]
    trans_vec = x[1:2]
    LazySets.translate!(body, trans_vec)

    return body
end

# # main
# atc_x = 0.25
# atc_y = 0.0
# bl = 0.75
# bw = 0.375

# # define body polygon at origin
# x0_min = atc_x - 1/2*bl
# x0_max = atc_x + 1/2*bl
# y0_min = atc_y - 1/2*bw
# y0_max = atc_y + 1/2*bw

# # ISSUE: will this be slow because it creates a new Polygon every time function is called?
# #   - probably, function has 32 allocations
# #   - rotate and translate functions might be allocating too, need to make in place
# #   - make sure function works properly first, then start optimizing

# #   - can define body0 as part of the veh struct
# body = VPolygon([[x0_min, y0_min], [x0_max, y0_min], [x0_max, y0_max], [x0_min, y0_max]])

# x = [1.0, 1.5, -1/8*pi]

# @btime body_rt = state_to_body(x, body)

# p1 = plot(body_rt, aspect_ratio=:equal)

# display(p1)