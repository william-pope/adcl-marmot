# Hamilton-Jacobi-Bellman demonstration

using Interpolations
using Plots

# notation:
#   x_p - horizontal position
#   y_p - vertical position
#   x - point [x_p, y_p]
#   y - pose [x, theta]
#   z - path, parameterized by time
#   sigma - speed control
#   a - steering control
#   T - target set
#   f - dynamics
#   rho - minimum turning radiuss
#   u - value function for given dynamics and target set
#   u_x, u_y, u_theta - partial derivatives of value function wrt x, y, theta

# TO-DO:
# - visualize somehow to verify (make 3D shape)
# - trace optimal paths using characteristic curves
# - make obstacle cost higher, just adjust color scale on plot

# 1) DEFINITIONS --- --- ---

# define discretized grid struct
struct StateGrid
    h_xy::Float64
    h_theta::Float64
    x_grid::StepRangeLen
    y_grid::StepRangeLen
    theta_grid::StepRangeLen
end

struct Vehicle
    c_vf::Float64   # max forward speed [m/s]
    c_vb::Float64   # max backward speed [m/s]
    c_phi::Float64  # max steering angle [rad]
    l::Float64      # wheelbase [m]
end

# finite difference approximations
# (eq 30)
function G_p1_ijk(theta_k, u_ip, u_jp, u_kp, u_kn, sg::StateGrid, v::Vehicle)
    rho = v.l/tan(v.c_phi)

    num = sg.h_xy/v.c_vf + abs(cos(theta_k))*u_ip + abs(sin(theta_k))*u_jp + sg.h_xy/(rho*sg.h_theta)*min(u_kp, u_kn)
    den = abs(cos(theta_k)) + abs(sin(theta_k)) + sg.h_xy/(rho*sg.h_theta)
    G_p1 = num/den

    return G_p1
end

# (eq 31)
function F_p1_ijk(theta_k, u_ip, u_jp, sg::StateGrid, v::Vehicle)
    num = sg.h_xy/v.c_vf + abs(cos(theta_k))*u_ip + abs(sin(theta_k))*u_jp
    den = abs(cos(theta_k)) + abs(sin(theta_k))
    F_p1 = num/den

    return F_p1
end

# (eq 32)
function G_n1_ijk(theta_k, u_in, u_jn, u_kp, u_kn, sg::StateGrid, v::Vehicle)
    rho = v.l/tan(v.c_phi)

    num = sg.h_xy/v.c_vb + abs(cos(theta_k))*u_in + abs(sin(theta_k))*u_jn + sg.h_xy/(rho*sg.h_theta)*min(u_kp, u_kn)
    den = abs(cos(theta_k)) + abs(sin(theta_k)) + sg.h_xy/(rho*sg.h_theta)
    G_n1 = num/den

    return G_n1
end

# (eq 33)
function F_n1_ijk(theta_k, u_in, u_jn, sg::StateGrid, v::Vehicle)
    num = sg.h_xy/v.c_vb + abs(cos(theta_k))*u_in + abs(sin(theta_k))*u_jn
    den = abs(cos(theta_k)) + abs(sin(theta_k))
    F_n1 = num/den

    return F_n1
end

# target set checker
function in_target_set(y, T_set)
    if y[1] < T_set[1][1] || y[1] > T_set[1][2]
        return false
    elseif y[2] < T_set[2][1] || y[2] > T_set[2][2]
        return false
    elseif y[3] < T_set[3][1] || y[3] > T_set[3][2]
        return false
    else
        return true
    end
end

# obstacle set checker
function in_obstacle_set(y, O_set)
    for On_set in O_set
        if y[1] < On_set[1][1] || y[1] > On_set[1][2]
            continue
        elseif y[2] < On_set[2][1] || y[2] > On_set[2][2]
            continue
        else
            return true
        end
    end

    return false
end

# initialize value approximations
function initialize_value_array(sg::StateGrid)
    Up = ones(Float64, length(sg.x_grid), length(sg.y_grid), length(sg.theta_grid))

    for i in 1:length(sg.x_grid)
        for j in 1:length(sg.y_grid)
            for k in 1:length(sg.theta_grid)
                x_i = sg.x_grid[i]
                y_j = sg.y_grid[j]
                theta_k = sg.theta_grid[k]

                y_ijk = [x_i, y_j, theta_k]

                if in_target_set(y_ijk, T_set) == true
                    Up[i,j,k] = 0
                else
                    Up[i,j,k] = 100
                end
            end
        end
    end

    U = deepcopy(Up)

    return U, Up
end

# computes value update for grid point ijk using HJB finite difference scheme
function update_value(U, i, j, k, sg::StateGrid, v::Vehicle)
    theta_k = sg.theta_grid[k]

    xi_k = Int(sign(cos(theta_k)))
    nu_k = Int(sign(sin(theta_k)))

    ip = i + xi_k
    in = i - xi_k
    jp = j + nu_k
    jn = j - nu_k
    k == length(sg.theta_grid) ? kp = 1 : kp = k + 1
    k == 1 ? kn = length(sg.theta_grid) : kn = k - 1

    u_ip = U[ip,j,k] 
    u_in = U[in,j,k]
    u_jp = U[i,jp,k]
    u_jn = U[i,jn,k]
    u_kp = U[i,j,kp]
    u_kn = U[i,j,kn]

    # compute u_ijk update
    G_p1 = G_p1_ijk(theta_k, u_ip, u_jp, u_kp, u_kn, sg, v)
    F_p1 = F_p1_ijk(theta_k, u_ip, u_jp, sg, v)
    G_n1 = G_n1_ijk(theta_k, u_in, u_jn, u_kp, u_kn, sg, v)
    F_n1 = F_n1_ijk(theta_k, u_in, u_jn, sg, v)

    up_ijk = min(G_p1, F_p1, G_n1, F_n1, U[i,j,k])

    return up_ijk
end

# main function to iteratively calculate value function using HJB
function solve_HJB_PDE(du_tol, max_reps, sg::StateGrid, v::Vehicle)
    N_x = length(sg.x_grid)
    N_y = length(sg.y_grid)
    N_theta = length(sg.theta_grid)

    # define ijk iterators for Gauss-Seidel sweeping scheme
    GS_sweeps = [[2:N_x-1, 2:N_y-1, 1:N_theta-1],     # FFF
                [2:N_x-1, 2:N_y-1, reverse(1:N_theta-1)],  # FFB
                [2:N_x-1, reverse(2:N_y-1), 1:N_theta-1],  # FBF
                [reverse(2:N_x-1), 2:N_y-1, 1:N_theta-1],  # BFF
                [2:N_x-1, reverse(2:N_y-1), reverse(1:N_theta-1)], # FBB
                [reverse(2:N_x-1), reverse(2:N_y-1), 1:N_theta-1], # BBF
                [reverse(2:N_x-1), 2:N_y-1, reverse(1:N_theta-1)], # BFB
                [reverse(2:N_x-1), reverse(2:N_y-1), reverse(1:N_theta-1)]]   # BBB

    # initialize U and Up
    (U, Up) = initialize_value_array(sg)

    # main function loop
    du_max = maximum(Up)
    rep = 0
    while du_max > du_tol && rep < max_reps
        rep += 1
        println("\nrep: ", rep)

        for sweep in GS_sweeps
            for i in sweep[1]
                for j in sweep[2]
                    for k in sweep[3]
                        x_i = sg.x_grid[i]
                        y_j = sg.y_grid[j]
                        theta_k = sg.theta_grid[k]

                        y_ijk = [x_i, y_j, theta_k]
                        
                        if in_obstacle_set(y_ijk, O_set) == false
                            Up[i,j,k] = update_value(U, i, j, k, sg, v)
                        end
                    end
                end
            end
            Up[:,:,end] = deepcopy(Up[:,:,1])
            Up[:,:,end] = Up[end:-1:1,:,end]

            # compare U and Up to check convergence
            dU = Up - U
            du_max = maximum(abs.(dU))
            println(du_max)

            U = deepcopy(Up)
        end
    end

    return U
end

# interpolate U to return a value for any point (x,y,theta)
function interp_value(y, U, sg::StateGrid)
    # adjust theta within bounds
    if y[3] > pi
        y[3] -= 2*pi
    elseif y[3] < -pi
        y[3] += 2*pi
    end

    # implements trilinear interpolation from Wikipedia
    if y[1] in sg.x_grid
        i_0 = indexin(y[1], sg.x_grid)[1]
    else
        i_0 = searchsortedfirst(sg.x_grid, y[1]) - 1
    end

    if y[2] in sg.y_grid
        j_0 = indexin(y[2], sg.y_grid)[1]
    else
        j_0 = searchsortedfirst(sg.y_grid, y[2]) - 1
    end

    if y[3] in sg.theta_grid
        k_0 = indexin(y[3], sg.theta_grid)[1]
    else
        k_0 = searchsortedfirst(sg.theta_grid, y[3]) - 1
    end

    i_1 = i_0 + 1
    j_1 = j_0 + 1
    k_0 == length(sg.theta_grid) ? k_1 = 1 : k_1 = k_0 + 1  # ISSUE

    x_0 = sg.x_grid[i_0]
    y_0 = sg.y_grid[j_0]
    theta_0 = sg.theta_grid[k_0]    # ISSUE: k_0 = 0 somehow

    x_1 = sg.x_grid[i_1]
    y_1 = sg.y_grid[j_1]
    theta_1 = sg.theta_grid[k_1]

    x_d = (y[1] - x_0)/(x_1 - x_0)
    y_d = (y[2] - y_0)/(y_1 - y_0)
    theta_d = (y[3] - theta_0)/(theta_1 - theta_0)

    u_000 = U[i_0, j_0, k_0]
    u_100 = U[i_1, j_0, k_0]
    u_010 = U[i_0, j_1, k_0]
    u_110 = U[i_1, j_1, k_0]
    u_001 = U[i_0, j_0, k_1]
    u_101 = U[i_1, j_0, k_1]
    u_011 = U[i_0, j_1, k_1]
    u_111 = U[i_1, j_1, k_1]

    u_00 = u_000*(1 - x_d) + u_100*x_d
    u_01 = u_001*(1 - x_d) + u_101*x_d
    u_10 = u_010*(1 - x_d) + u_110*x_d
    u_11 = u_011*(1 - x_d) + u_111*x_d

    u_0 = u_00*(1 - y_d) + u_10*y_d
    u_1 = u_01*(1 - y_d) + u_11*y_d

    u_itp = u_0*(1 - theta_d) + u_1*theta_d

    return u_itp
end

function find_path(y_0, U, T_set, dt, sg::StateGrid, v::Vehicle)
    max_steps = 5000

    if typeof(y_0) != Vector{Float64}
        y_0 = convert(Vector{Float64}, y_0)
    end

    y_k = y_0
    y_path = [y_k]  

    u_path = []

    step = 0
    
    while in_target_set(y_k, T_set) == false && step < max_steps
        step += 1

        # calculate optimal action
        u_k = optimal_action(y_k, U, sg, v)
        push!(u_path, u_k)

        # simulate forward one time step
        y_k1 = runge_kutta_4(car_EoM, y_k, u_k, dt, v)
        push!(y_path, y_k1)

        y_k = deepcopy(y_k1)
    end

    println("steps in path: ", step)

    return y_path, u_path, step
end

# ISSUE: problem with steering (dt too large?, interp not good?)
#   - creates repeating half-circle pattern, keeps steering past correct direction
#   - multiple time steps in each half-circle, dt seems fine
#   - value function seems fine

#   - could be an issue with steering around u_theta = 0
#   - create separate function for optimal action solver
#   - du_dtheta doesn't like to settle to 0 (it should), might need finer theta grid

#   - (***) issue navigating through saddles/minima in value function

#   - two notable behaviors:
#       - car keeps steering past what should be the optimal direction (stuck in turn)
#       - car stops/reverses at +1/-1 speed saddle (equal value both directions)

#   - works great with no obstacles (because no saddles/minima)

#   - sorta working with a symetric obstacle (?)

#   - saddles should be unstable, but weird steering mechanics may change that
#   - there are some weird minima behind obstacles (!)
#   - ISSUE WITH VALUE FUNCTION AND OBSTACLES

#   - assymetric obstacles:
#       - U gets stuck temporarily at certain values during computation

#   - when max_steps reached for seemingly short path, it's because speed is switching
#       back and forth, while steering remains fixed to one side

#   LOOK AT FURTHER:
#   - speed input chooses first, based on value gradient at the given theta, and always chooses
#       speed direction that will push downhill *within the current theta slice*
#   - may cause issues where best action in 3D gradient is good in theta-dimension, but bad in
#       xy-dimension. This action won't be chosen, because speed input chooses before steering input

#   - may have separate/similar issue with chattering around saddle/minima

#   - issue with theta wrapping when calculating action/du_dtheta?
#       - seems like it has trouble when theta=+pi
#       - RK4 simulates to theta outside bounds, then doesn't know what to do?
#       - this solved some problems, but still issues around +/-pi in certain areas of workspace
#       - seems to occur when approaching +pi from <+pi, and driving in reverse

#   - did changing fwd/bkwd speed get rid of minima behind obstacle? not a great solution if so...

# NEED TO UNDERSTAND FORWARD/BACKWARD CHATTERING PROBLEM

# calculates optimal action as a function of state
function optimal_action(y, U, sg::StateGrid, v::Vehicle)
    # adjust position when near edge for boundary problems
    if abs(y[1] - sg.x_grid[1]) < 2*sg.h_xy
        y[1] = sg.x_grid[3]
    elseif abs(y[1] - sg.x_grid[end]) < 2*sg.h_xy
        y[1] = sg.x_grid[end-2]
    end

    if abs(y[2] - sg.y_grid[1]) < 2*sg.h_xy
        y[2] = sg.y_grid[3]
    elseif abs(y[2] - sg.y_grid[end]) < 2*sg.h_xy
        y[2] = sg.y_grid[end-2]
    end

    # approximate partial derivatives using centered difference scheme
    du_dx = (interp_value(y+[sg.h_xy,0,0], U, sg) - interp_value(y-[sg.h_xy,0,0], U, sg))/(2*sg.h_xy)
    du_dy = (interp_value(y+[0,sg.h_xy,0], U, sg) - interp_value(y-[0,sg.h_xy,0], U, sg))/(2*sg.h_xy)
    du_dtheta = (interp_value(y+[0,0,sg.h_theta], U, sg) - interp_value(y-[0,0,sg.h_theta], U, sg))/(2*sg.h_theta)

    # calculate optimal action
    A = [[v.c_vf, 0],
        [v.c_vf, v.c_phi],
        [v.c_vf, -v.c_phi],
        [-v.c_vb, 0],
        [-v.c_vb, v.c_phi],
        [-v.c_vb, -v.c_phi]]

    HJB_min(a) = a[1]*(du_dx*cos(y[3]) + du_dy*sin(y[3]) + du_dtheta*(1/v.l)*tan(a[2]))
    a_opt = argmin(HJB_min, A)

    return a_opt
end

# equations of motion for 3 DoF kinematic bicycle model
function car_EoM(x, u, param)
    x_dot = [u[1]*cos(x[3]),
            u[1]*sin(x[3]),
            u[1]*(1/param.l)*tan(u[2])]

    return x_dot
end

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(EoM::Function, x_k, u, dt, param)
    w1 = EoM(x_k, u, param)
    w2 = EoM(x_k + w1*dt/2, u, param)
    w3 = EoM(x_k + w2*dt/2, u, param)
    w4 = EoM(x_k + w3*dt, u, param)

    x_k1 = x_k + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    # adjust theta within bounds
    if x_k1[3] > pi
        x_k1[3] -= 2*pi
    elseif x_k1[3] < -pi
        x_k1[3] += 2*pi
    end

    return x_k1
end


# 2) PARAMETERS --- --- ---
# parameters that determine behavior:
#   - h_xy
#   - h_theta
#   - dt

# vehicle parameters
marmot = Vehicle(1.5, 0.75, 0.475, 0.324)   
unit_car = Vehicle(1.0, 0.95, 0.5, 1.0)   

# define workspace
W_x = [-10, 10]
W_y = [-10, 10]
W_theta = [-pi, pi]
W_set = [W_x, W_y, W_theta]

# define target region
T_x = [-1, 1]
T_y = [4, 6]
T_theta = [-pi, pi]
T_set = [T_x, T_y, T_theta]

# define obstacles
O1_x = [-6, 2]
O1_y = [-3, -2]
O1_set = [O1_x, O1_y]

O_set = [O1_set]
# O_set = []

# initialize state grid /// /// /// /// /// /// ///
my_h_xy = 0.5
my_h_theta = deg2rad(5)

sg = StateGrid(my_h_xy, 
                my_h_theta,
                W_x[1]:my_h_xy:W_x[2],
                W_y[1]:my_h_xy:W_y[2],
                W_theta[1]:my_h_theta:W_theta[2])


# 3) MAIN --- --- ---
println("\nstart --- --- ---")

du_tol = 0.01
max_reps = 100
@time U = solve_HJB_PDE(du_tol, max_reps, sg, unit_car)


# 4) PLOTS --- --- ---
W_x_pts = [W_x[1],W_x[2],W_x[2],W_x[1],W_x[1]]
W_y_pts = [W_y[1],W_y[1],W_y[2],W_y[2],W_y[1]]

T_x_pts = [T_x[1],T_x[2],T_x[2],T_x[1],T_x[1]]
T_y_pts = [T_y[1],T_y[1],T_y[2],T_y[2],T_y[1]]

O1_x_pts = [O1_x[1],O1_x[2],O1_x[2],O1_x[1],O1_x[1]]
O1_y_pts = [O1_y[1],O1_y[1],O1_y[2],O1_y[2],O1_y[1]]

# plot u(x,y,0) as color map
p1_k = Int(round(length(sg.theta_grid)/2,digits=0))
p1 = heatmap(sg.x_grid, sg.y_grid, transpose(U[:,:,p1_k]), clim=(0,20),
            aspect_ratio=:equal, size=(700,600),
            xlabel="x-axis", ylabel="y-axis", title="HJB Value Function: u(x, y, theta=0)")

plot!(p1, W_x_pts, W_y_pts, linewidth=2, linecolor=:black, label="Workspace")
plot!(p1, T_x_pts, T_y_pts, linewidth=2, linecolor=:green, label="Target Set")
display(p1)

# plot u(x,y,pi/2) as color map
p2_k = Int(round(length(sg.theta_grid)*3/4,digits=0))
p2 = heatmap(sg.x_grid, sg.y_grid, transpose(U[:,:,p2_k]), clim=(0,20),
            aspect_ratio=:equal, size=(700,600),
            xlabel="x-axis", ylabel="y-axis", title="HJB Value Function: u(x, y, theta=pi/2)")

plot!(p2, W_x_pts, W_y_pts, linewidth=2, linecolor=:black, label="Workspace")
plot!(p2, T_x_pts, T_y_pts, linewidth=2, linecolor=:green, label="Target Set")
display(p2)

# plot u(x,y,pi/2) as color map
p3_k = Int(round(length(sg.theta_grid)*5/8,digits=0))
p3 = heatmap(sg.x_grid, sg.y_grid, transpose(U[:,:,p3_k]), clim=(0,20),
            aspect_ratio=:equal, size=(700,600),
            xlabel="x-axis", ylabel="y-axis", title="HJB Value Function: u(x, y, theta=pi/4)")

plot!(p3, W_x_pts, W_y_pts, linewidth=2, linecolor=:black, label="Workspace")
plot!(p3, T_x_pts, T_y_pts, linewidth=2, linecolor=:green, label="Target Set")
display(p3)

# plot optimal path from y_0 to target set
p_path = plot(aspect_ratio=:equal, size=(700,600))
plot!(p_path, W_x_pts, W_y_pts, linewidth=2, linecolor=:black, label="Workspace")
plot!(p_path, T_x_pts, T_y_pts, linewidth=2, linecolor=:green, label="Target Set")
plot!(p_path, O1_x_pts, O1_y_pts, linewidth=2, linecolor=:red, label="Obstacle 1")

dt = 0.1

# wrap_around = [[-8, -8, pi/2],
#                 [-4, -8, pi/2],
#                 [0, -8, pi/2],
#                 [4, -8, pi/2],
#                 [8, -8, pi/2],
#                 [8, -4, pi/2],
#                 [8, 0, pi/2],
#                 [8, 4, pi/2],
#                 [8, 8, pi/2],
#                 [4, 8, pi/2],
#                 [0, 8, pi/2],
#                 [-2, 8, pi/2],
#                 [-4, 8, pi/2],
#                 [-6, 8, pi/2],
#                 [-8, 8, pi/2],
#                 [-8, 6, pi/2],
#                 [-8, 4, pi/2],
#                 [-8, 2, pi/2],
#                 [-8, 0, pi/2],
#                 [-8, -4, pi/2]]

asym_obs = [[-8, -8, pi/2],
            [-6, -8, pi/2],
            [-4, -8, pi/2],
            [-2, -8, pi/2],
            [0, -8, pi/2],
            [2, -8, pi/2],
            [4, -8, pi/2],
            [-5, -6, pi/2],
            [-3, -6, pi/2],
            [-1, -6, pi/2],
            [1, -6, pi/2],
            [-4, -4.5, pi/2],
            [-2, -4.5, pi/2],
            [0, -4.5, pi/2]]

initial_poses = asym_obs

for y_0 in initial_poses
    (y_path, u_path) = find_path(y_0, U, T_set, dt, sg, unit_car)
    plot!(p_path, getindex.(y_path,1), getindex.(y_path,2), linewidth=2, label="")
end

# for _ in 1:50
#     y_0 = [rand(W_x[1]:0.01:W_x[2]), rand(W_y[1]:0.01:W_y[2]), rand(-pi:0.01:pi)]
#     (y_path, u_path, steps) = find_path(y_0, U, T_set, dt, sg, unit_car)
#     if steps == 2000
#         println(y_path[end])
#         println(u_path[end])
#         println("")
#         plot!(p_path, getindex.(y_path,1), getindex.(y_path,2), linewidth=2, color=:red, label="")
#     else
#         plot!(p_path, getindex.(y_path,1), getindex.(y_path,2), linewidth=2, color=:blue, label="")
#     end
# end

display(p_path)
