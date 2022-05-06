# final project

using Plots
using JLD

include("HJB_functions.jl")

# need for MDP:
#   - state space (S)
#   - action space (A)
#   - transition function (T(sp|s,a))
#   - reward function (R(s,a))
#   - discount (gamma)

# TO-dO:
#   - deal with a mountain of bugs lol

function continuous_MCTS(s_0, V_HJB, dt, S, A, gamma, std_v, std_phi, sims, w_max, d_max, c_UCB, EoM::Function, sg::StateGrid, veh::Vehicle)
    s_type = typeof(S[1])
    a_type = typeof(A[1])

    N = Dict{Tuple{s_type, a_type}, Int}()
    S_p = Dict{Tuple{s_type, a_type}, Vector{s_type}}()
    Q = Dict{Tuple{s_type, a_type}, Float64}()

    for a in A
        N[(s_0, a)] = 0
        S_p[(s_0, a)] = []
        Q[(s_0, a)] = 0.0 
    end 

    # conduct MCTS to fill Q_sa dictionary
    for sim in 1:sims
        s = s_0
        s_ro = s_0
        d = 1

        # SEARCH ---
        at_edge = false
        while at_edge == false
            # println("s: ", s)

            # if max depth is reached
            if d >= d_max 
                s_ro = deepcopy(s)     
                at_edge = true
            end 

            # if not at max depth, keep stepping through tree until new node found
            a_UCB = select_UCB_action(s, c_UCB, A, N, Q)
            # println("a_UCB: ", a_UCB)

            # if max width not yet reached for given layer
            if size(S_p[(s, a_UCB)], 1) < w_max
                # println("layer not filled, propagating new state node")
                s_p = gen_state(EoM, s, a_UCB, std_v, std_phi, dt, veh)
                # println("s_p: ", s_p)
                
                push!(S_p[(s, a_UCB)], s_p)

                for a in A
                    N[(s_p, a)] = 0
                    S_p[(s_p, a)] = []
                    Q[(s_p, a)] = 0.0
                end 

                # println("S_p dict: ")
                # display(S_p)

                s_ro = deepcopy(s_p)
                at_edge = true
            else
                # println("layer filled, choosing node and continuing to step")
                s_p = rand(S_p[(s, a_UCB)])
                s = deepcopy(s_p)

                if d > 2
                    println(d)
                end

                d += 1
            end
        end

        # ROLLOUT ---
        v_ro_est = interp_value(s_ro, V_HJB, sg)
        # println("v_ro_est from s_p: ", v_ro_est)

        # BACKUP ---
        q_sa_p = v_ro_est
        s_p = deepcopy(s_ro)
        # println("\nbackup:")

        at_root = false
        while at_root == false
            # println("s_p: ", s_p)
            key_sa = [key for (key, val) in S_p if s_p in val]  
            
            s = key_sa[1][1]
            a_UCB = key_sa[1][2]
            
            # println("a_UCB: ", a_UCB)
            # println("s: ", s)

            r_sa = gen_reward(s, a_UCB, dt, T_xy_set, T_theta_set, veh)
            q_sa = r_sa + gamma*q_sa_p

            N[(s, a_UCB)] += 1
            Q[(s, a_UCB)] += (q_sa - Q[(s, a_UCB)])/N[(s, a_UCB)]

            if s == s_0
                # println("backed up to root")
                at_root = true
            end

            s_p = deepcopy(s)
            q_sa_p = deepcopy(q_sa)
        end

        # println("Q dict: ")
        # display(Q)
    end 

    # chooses best action for current root node
    a_best = A[1]
    Q_best = -Inf
    for a in A
        Q_sa = Q[(s_0, a)]
        if Q_sa > Q_best
            a_best = a
            Q_best = Q_sa
        end
    end

    # display(Q)
    # println(length(N))

    return a_best, N, S_p
end 

function select_UCB_action(s, c_UCB, A, N, Q) 
    N_s = sum(N[(s,a)] for a in A)

    a_UCB = A[1]
    UCB_max = -Inf
    for a in A
        N_sa = N[(s,a)]

        if N_sa == 0
            UCB = Inf
        else
            exp_term = c_UCB*sqrt(log(N_s)/N_sa)
            UCB = Q[(s,a)] + exp_term
        end

        if UCB > UCB_max
            a_UCB = a
            UCB_max = UCB
        end
    end

    return a_UCB
end 

function gen_state(car_EoM::Function, s, a_UCB, std_v, std_phi, dt, veh::Vehicle)
    # simulate forward one time step
    s_p = runge_kutta_4(car_EoM, s, a_UCB, std_v, std_phi, dt, veh)

    return s_p
end

function gen_reward(s, a, dt, T_xy_set, T_theta_set, veh::Vehicle)
    if in_target_set(s, T_xy_set, T_theta_set, veh) == true
        r_sa = 0
    else
        r_sa = -dt
    end

    return r_sa
end

function mcts_planner(s_0, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, EoM::Function, sg::StateGrid, veh::Vehicle)
    # ISSUE: does number of sims need to be direct multiple of size of A?
    d_max = 10
    c_UCB = 8

    max_steps = 5000

    a_hist_mcts = []
    s_hist_mcts = [s_0]

    p = true

    s = s_0
    step = 0
    while in_target_set(s, T_xy_set, T_theta_set, veh) == false && step < max_steps
        step += 1

        a_mcts, N, S_p = continuous_MCTS(s, V_HJB, dt, S, A, gamma, std_v, std_phi, sims, w_max, d_max, c_UCB, EoM, sg, veh)

        s_p = runge_kutta_4(car_EoM, s, a_mcts, std_v, std_phi, dt, unit_car)

        # println("step, s_p: ", (step, s_p))

        # if p == true
        #     display(S_p[([-4,-8,pi/2],[1.0,0.0])])
        #     p = false
        # end

        push!(a_hist_mcts, a_mcts)
        push!(s_hist_mcts, s_p)

        s = deepcopy(s_p)
    end

    println("steps in MCTS path: ", step)

    return s_hist_mcts, a_hist_mcts, step
end


# 2) PARAMETERS --- --- ---

# vehicle parameters
unit_car = Vehicle(1.0, 0.75, 0.5, 0.5, 0.75, 0.25, 0.125)   

# define workspace
W_set = [[-10.0 -10.0];
        [10.0 -10.0];
        [10.0 10.0];
        [-10.0 10.0]]

# define target set
T_xy_set = [[-1.0 4.0];
            [1.0 4.0];
            [1.0 6.0];
            [-1.0 6.0]]

T_theta_set = [[-pi, pi]]

# define obstacles
O1_set = [[-2.0 -1.0];
        [3.5 -1.0];
        [3.5 1.0];
        [-2.0 1.0]]

O2_set = [[-6.5 -5.0];
        [-2.0 -5.0];
        [-2.0 -3.0];
        [-6.5 -3.0]]

O_set = [O1_set, O2_set]
# O_set = []

# initialize state grid
h_xy = 0.125
h_theta = deg2rad(2.5)

sg = StateGrid(h_xy, 
                h_theta,
                minimum(W_set[:,1]) : h_xy : maximum(W_set[:,1]),
                minimum(W_set[:,2]) : h_xy : maximum(W_set[:,2]),
                -pi : h_theta : pi)


# 3) MAIN --- --- ---
println("\n--- START ---")

# HJB
run_HJB = true

if run_HJB == true
    du_tol = 0.01
    max_reps = 100
    @time U = solve_HJB_PDE(du_tol, max_reps, sg, O_set, unit_car)
    V_HJB = -U

    n_nodes = size(sg.x_grid, 1) * size(sg.y_grid, 1) * size(sg.theta_grid, 1)
    println("total number of grid nodes: ", n_nodes)

    save("V_HJB.jld", "V_HJB", V_HJB)
else
    V_HJB = load("V_HJB.jld", "V_HJB")
    U = -V_HJB
end

# MCTS
S = [[-10.0, 10.0],
    [-10.0, 10.0],
    [-pi, pi]]

A_v = [-0.751, 1.0]
A_phi = [-0.5, 0.0, 0.5]
A = vec([[a_v, a_phi] for a_phi in A_phi, a_v in A_v])

gamma = 0.95

# ISSUE: noise works, but should add to action instead of x_dot
#   - needs different scaling for x,y (orientation) and theta (units)
#   - should be a "correct" way to add noise to bicycle model

# not doing anything?

s_0_mcts = [-4, -8, -0.55*pi]
dt = 0.1

println("planning path 1")
std_v = 0.02
std_phi = 0.01
sims = 324
w_max = 3

s_0_mcts = [-5, -7, -0.5*pi]

s_hist1_mcts, a_hist1_mcts = mcts_planner(s_0_mcts, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, car_EoM, sg, unit_car)

println("planning path 2")
std_v = 0.1
std_phi = 0.05
sims = 324
w_max = 3

s_0_mcts = [-3.75, -8.25, -0.5*pi]

s_hist2_mcts, a_hist2_mcts = mcts_planner(s_0_mcts, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, car_EoM, sg, unit_car)

println("planning path 3")
std_v = 0.5
std_phi = 0.25
sims = 324
w_max = 3

s_0_mcts = [-2.5, -9.5, -0.5*pi]

s_hist3_mcts, a_hist3_mcts = mcts_planner(s_0_mcts, V_HJB, dt, S, A, gamma, std_v, std_phi, w_max, sims, T_xy_set, T_theta_set, car_EoM, sg, unit_car)

# s_0_hjb = [-4.0, -8, pi/2]
# y_path_HJB, u_path_HJB = HJB_planner(s_0_hjb, U, T_xy_set, T_theta_set, dt, sg, unit_car)

# display(s_hist_mcts)
# display(a_hist_mcts)


# 4) PLOTS --- --- ---

# plot U as heat map
if run_HJB == true
    for k_plot in LinRange(1, size(sg.theta_grid, 1), 15)
        k_plot = Int(round(k_plot, digits=0))
        theta_k = round(rad2deg(sg.theta_grid[k_plot]), digits=3)

        p_k = heatmap(sg.x_grid, sg.y_grid, transpose(V_HJB[:,:,k_plot]), clim=(-20,0),
                    aspect_ratio=:equal, size=(675,600),
                    xlabel="x-axis [m]", ylabel="y-axis [m]",
                    right_margin = 4Plots.mm,
                    top_margin = -8Plots.mm,
                    bottom_margin = -8Plots.mm)
                    # , title="HJB Value Function: u(x, y, theta=$theta_k)"

        plot_polygon(W_set, p_k, 3, :black, "Workspace")
        plot_polygon(T_xy_set, p_k, 3, :green, "Target Set")
        plot_polygon(O1_set, p_k, 3, :red, "Obstacle")
        plot_polygon(O2_set, p_k, 3, :red, "")

        display(p_k)
    end
end

# plot optimal path from y_0 to target set
p_path_mcts = plot(aspect_ratio=:equal, size=(600,600), legend=:topright)
            # xlabel="x-axis [m]", ylabel="y-axis [m]")

plot_polygon(W_set, p_path_mcts, 3, :black, "Workspace")
plot_polygon(T_xy_set, p_path_mcts, 3, :green, "Target Set")
plot_polygon(O1_set, p_path_mcts, 3, :red, "Obstacle")
plot_polygon(O2_set, p_path_mcts, 3, :red, "")

# plot!(p_path_mcts, getindex.(y_path_HJB,1), getindex.(y_path_HJB,2), 
#     color=:black, linewidth=2.5, linestyle=:dash, label="HJB Path")

plot!(p_path_mcts, getindex.(s_hist1_mcts,1), getindex.(s_hist1_mcts,2),
    color=:green, linewidth=2.5, label="MCTS, σ1")

plot!(p_path_mcts, getindex.(s_hist2_mcts,1), getindex.(s_hist2_mcts,2),
    color=:blue, linewidth=2.5, label="MCTS, σ2")

plot!(p_path_mcts, getindex.(s_hist3_mcts,1), getindex.(s_hist3_mcts,2),
    color=:purple, linewidth=2.5, label="MCTS, σ3")

display(p_path_mcts)