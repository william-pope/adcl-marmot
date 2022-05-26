# MCTS functions

using Plots

algs_path_mac = "/Users/willpope/Desktop/Research/marmot-algs/"
include(algs_path_mac*"HJB-planner/HJB_functions.jl")

struct Solver
    sims::Int
    w_max::Int
    d_max::Int
    c_UCB::Float64
end

function continuous_MCTS(s_0, V_HJB, Dt, K_sub, A, gamma, std_v, std_phi, slv::Solver, EoM::Function, env::Environment, veh::Vehicle)
    s_type = typeof(s_0)
    a_type = typeof(A[1])

    S = Dict{Int, Tuple{s_type, Int}}()
    IA = Dict{Int, Tuple{Int, a_type}}()
    I_p = Dict{Tuple{Int, a_type}, Vector{Int}}()
    Subpaths = Dict{Tuple{Int, Int}, Vector{s_type}}()
    N = Dict{Tuple{Int, a_type}, Int}()
    Q = Dict{Tuple{Int, a_type}, Float64}()
    R = Dict{Tuple{Int, a_type, Int}, Float64}()

    N_hist = Dict{Tuple{Int, a_type}, Vector{Int}}()
    Q_hist = Dict{Tuple{Int, a_type}, Vector{Float64}}()
    for a in A
        N_hist[(1, a)] = []
        Q_hist[(1, a)] = []
    end

    S[1] = (s_0, 0)
    for a in A
        I_p[(1, a)] = []
        N[(1, a)] = 0
        Q[(1, a)] = 0.0     # ?: does this impact anything? since 0.0 is actually higher than most
    end 

    # p_nodes = plot(aspect_ratio=:equal, size=(800,600))
    # plot!(p_nodes, [s_0[1]], [s_0[2]], 
    #     markersize=6, markershape=:circle, markerstrokewidth=1, markercolor=:black,
    #     label="")

    # println("s_0: ", s_0)

    # conduct MCTS to fill Q_sa dictionary
    for sim in 1:slv.sims
        # println("sim: ", sim)
        
        i = 1
        i_ro = 1
        d = 1

        # SEARCH ---
        at_edge = false
        while at_edge == false
            # println("i: ", i)

            # if max depth is reached
            if d >= slv.d_max 
                i_ro = deepcopy(i)     
                at_edge = true
            end 

            # if not at max depth, keep stepping through tree until new node found
            a_UCB = select_UCB_action(i, slv.c_UCB, A, N, Q)
            # println("a_UCB: ", a_UCB)

            # if max width not yet reached for given layer
            if size(I_p[(i, a_UCB)], 1) < slv.w_max
                # layer not filled, propagating new state node
                s_p, subpath_p = gen_state(S[i][1], a_UCB, std_v, std_phi, Dt, K_sub, EoM, veh)
                i_p = maximum(keys(S)) + 1

                S[i_p] = (s_p, d)
                IA[i_p] = (i, a_UCB)
                push!(I_p[(i, a_UCB)], i_p) 
                Subpaths[i, i_p] = subpath_p

                r_subpath = gen_path_reward(subpath_p, a_UCB, Dt, K_sub, env, veh)
                R[(i, a_UCB, i_p)] = r_subpath

                # println("s_p: ", s_p)
                # println("r_subpath: ", r_subpath)
                # println("")

                if a_UCB[1] > 0
                    mcolor = :blue
                else
                    mcolor = :red
                end

                # plot!(p_nodes, [s_p[1]], [s_p[2]], 
                #     markersize=3, markershape=:circle, markercolor=mcolor, 
                #     markerstrokewidth=1,
                #     label="")
                
                # plot!(p_nodes, [getindex.(subpath_p,1); s_p[1]], [getindex.(subpath_p,2); s_p[2]],
                #     linecolor=mcolor,
                #     label="")

                # display(p_nodes)

                for a in A
                    I_p[(i_p, a)] = []
                    N[(i_p, a)] = 0
                    Q[(i_p, a)] = 0.0
                end 

                i_ro = deepcopy(i_p)
                at_edge = true
            else
                # layer filled, choosing node with a_UCB and continuing to step
                i_p = rand(I_p[(i, a_UCB)])
                i = deepcopy(i_p)

                d += 1
            end
        end

        # ROLLOUT ---
        v_ro_est = interp_value(S[i_ro][1], V_HJB, env) 

        # BACKUP ---
        q_sa_p = v_ro_est
        i_p = deepcopy(i_ro)

        at_root = false
        while at_root == false
            (i, a_UCB) = IA[i_p]

            r_sa = R[(i, a_UCB, i_p)] 
            q_sa = r_sa + gamma*q_sa_p
            # ?: shouldn't this be using V(s'), not Q(s',a_UCB)?
            #   - Q(s,a) = R(s,a) + γ*E[V(s')]
            #   - V(s') = future reward when following some policy (optimal?)

            N[(i, a_UCB)] += 1
            Q[(i, a_UCB)] += (q_sa - Q[(i, a_UCB)])/N[(i, a_UCB)]

            if i == 1
                # backed up to root
                at_root = true
            end

            i_p = deepcopy(i)
            q_sa_p = deepcopy(q_sa)
        end

        for a in A
            push!(N_hist[(1, a)], N[(1, a)])
            push!(Q_hist[(1, a)], Q[(1, a)])
        end
    end 

    p_Q_hist = plot(xlabel="sim", ylabel="Q(s0,a)", legend=:bottomleft)
    for ai in 1:size(A,1)
        plot!(p_Q_hist, 1:slv.sims, getindex.(Q_hist[(1, A[ai])]),
            label="a$ai")
    end
    display(p_Q_hist)

    p_N_hist = plot(xlabel="sim", ylabel="N(s0,a)", legend=:topleft)
    for ai in 1:size(A,1)
        plot!(p_N_hist, 1:slv.sims, getindex.(N_hist[(1, A[ai])]),
            label="a$ai")
    end
    display(p_N_hist)
    
    # chooses best action for current root node
    a_best = A[1]
    Q_best = -Inf
    for a in A
        Q_sa = Q[(1, a)]
        if Q_sa > Q_best
            a_best = a
            Q_best = Q_sa
        end
    end

    # i_max = maximum(keys(S))
    # println("\nnumber of children at each node:")
    # for i in 1:i_max
    #     ip_ct = 0
    #     for a in A
    #         for i_p in I_p[(i,a)]
    #             if isempty(i_p) == false
    #                 ip_ct += 1
    #             end
    #         end
    #     end
    #     println("i = ", i, ",  ", "\td = ", S[i][2], ", \tc = ", ip_ct)
    # end
    # println("")

    for a in A
        println("a: ", a, ",  q(s0,a) = ", Q[(1, a)])
    end

    println("a_best: ", a_best)
    println("")

    return a_best, S, Subpaths
end 

function select_UCB_action(i, c_UCB, A, N, Q) 
    N_i = sum(N[(i,a)] for a in A)

    a_UCB = A[1]
    UCB_max = -Inf
    for a in shuffle(A)
        N_ia = N[(i,a)]

        if N_ia == 0
            UCB = Inf
        else
            exp_term = c_UCB*sqrt(log(N_i)/N_ia)
            UCB = Q[(i,a)] + exp_term
        end

        if UCB > UCB_max
            a_UCB = a
            UCB_max = UCB
        end
    end

    return a_UCB
end 

function gen_state(s, a_UCB, std_v, std_phi, Dt, K_sub, car_EoM::Function, veh::Vehicle)
    dt = Dt/K_sub

    subpath_p = []   # sub_path should not contain final node
    s_k = s
    s_k1 = s
    for _ in 1:K_sub
        push!(subpath_p, s_k)

        s_k1 = runge_kutta_4(s_k, a_UCB, std_v, std_phi, dt, car_EoM, veh)
        s_k = deepcopy(s_k1)
    end

    s_p = s_k1

    return s_p, subpath_p
end

function gen_path_reward(subpath_p, a, Dt, K_sub, env::Environment, veh::Vehicle)
    dt = Dt/K_sub
    r_path = 0

    for s_sub in subpath_p
        # event rewards
        if in_obstacle_set(s_sub, env, veh) == true || in_workspace(s_sub, env, veh) == false
            r_path += -50.0 
            break   # doesn't evaluate path after collision
        end

        # duration rewards
        if in_target_set(s_sub, env, veh) == false
            r_path += -dt
        end
    end

    return r_path
end

function mcts_planner(s_0, V_HJB, Dt, K_sub, S, A, gamma, std_v, std_phi, slv::Solver, EoM::Function, env::Environment, veh::Vehicle)
    max_steps = 5000

    a_hist_mcts = []
    s_hist_mcts = [s_0]

    s = s_0
    step = 0
    while in_target_set(s, env, veh) == false && step < max_steps
        step += 1

        a_mcts = continuous_MCTS(s, V_HJB, Dt, K_sub, S, A, gamma, std_v, std_phi, slv, EoM, env, veh)

        s_p = runge_kutta_4(s, a_mcts, std_v, std_phi, Dt, K_sub, car_EoM, unit_car)   # needs K_sub

        push!(a_hist_mcts, a_mcts)
        push!(s_hist_mcts, s_p)

        s = deepcopy(s_p)
    end

    println("steps in MCTS path: ", step)

    return s_hist_mcts, a_hist_mcts, step
end



#   - need to figure out what to do when path violates
#       - theory of sampling-based MDP: collect values from possible outcomes
#       - so if invalid, should return some bad value up tree, which will be weighted by number of times it occurs in sampling

#       - should probably truncate branch after obstacle reward has been backed up, don't want to continue search
#           - UCB should handle this, but sorta leaves it to be determined by reward structure
#           - ?: think it's fine to insert a manual hard constraint on this, although may need some more thinking (?)
#           - could also return that sp is N/A, or first point in path that collides (not sure if matters)

#       - check Himanshu thesis to align approaches