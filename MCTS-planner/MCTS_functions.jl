# MCTS functions

struct Solver
    sims::Int
    w_max::Int
    d_max::Int
    c_UCB::Float64
end

function continuous_MCTS(s_0, V_HJB, dt, S, A, gamma, std_v, std_phi, slv::Solver, EoM::Function, env::Environment, veh::Vehicle)
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
    for sim in 1:slv.sims
        s = s_0
        s_ro = s_0
        d = 1

        # SEARCH ---
        at_edge = false
        while at_edge == false
            # println("s: ", s)

            # if max depth is reached
            if d >= slv.d_max 
                s_ro = deepcopy(s)     
                at_edge = true
            end 

            # if not at max depth, keep stepping through tree until new node found
            a_UCB = select_UCB_action(s, slv.c_UCB, A, N, Q)
            # println("a_UCB: ", a_UCB)

            # if max width not yet reached for given layer
            if size(S_p[(s, a_UCB)], 1) < slv.w_max
                # println("layer not filled, propagating new state node")
                s_p = gen_state(s, a_UCB, std_v, std_phi, dt, EoM, veh)
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

                d += 1
            end
        end

        # ROLLOUT ---
        v_ro_est = interp_value(s_ro, V_HJB, env)   # using HJB value function
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

            r_sa = gen_reward(s, a_UCB, dt, env, veh)
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
    a_best = rand(A)
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

    return a_best
end 

function select_UCB_action(s, c_UCB, A, N, Q) 
    N_s = sum(N[(s,a)] for a in A)

    a_UCB = rand(A)
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

function gen_state(s, a_UCB, std_v, std_phi, dt, car_EoM::Function, veh::Vehicle)
    # simulate forward one time step
    s_p = runge_kutta_4(car_EoM, s, a_UCB, std_v, std_phi, dt, veh)

    return s_p
end

function gen_reward(s, a, dt, env::Environment, veh::Vehicle)
    if in_target_set(s, env, veh) == true
        r_sa = 0
    else
        r_sa = -dt
    end

    return r_sa
end

function mcts_planner(s_0, V_HJB, dt, S, A, gamma, std_v, std_phi, slv::Solver, EoM::Function, env::Environment, veh::Vehicle)
    # Q: does number of sims need to be direct multiple of size of A?
    max_steps = 5000

    a_hist_mcts = []
    s_hist_mcts = [s_0]

    s = s_0
    step = 0
    while in_target_set(s, env, veh) == false && step < max_steps
        step += 1

        a_mcts = continuous_MCTS(s, V_HJB, dt, S, A, gamma, std_v, std_phi, slv, EoM, env, veh)

        s_p = runge_kutta_4(car_EoM, s, a_mcts, std_v, std_phi, dt, unit_car)

        push!(a_hist_mcts, a_mcts)
        push!(s_hist_mcts, s_p)

        s = deepcopy(s_p)
    end

    println("steps in MCTS path: ", step)

    return s_hist_mcts, a_hist_mcts, step
end