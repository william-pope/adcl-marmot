# HJB_utils.jl

using LazySets

include("dynamics_models.jl")

# 4th-order Runge-Kutta integration scheme
function runge_kutta_4(x_k::Vector{Float64}, u_k::Vector{Float64}, dt, EoM, veh, sg)
    w1 = EoM(x_k, u_k, veh)
    w2 = EoM(x_k + w1*dt/2, u_k, veh)
    w3 = EoM(x_k + w2*dt/2, u_k, veh)
    w4 = EoM(x_k + w3*dt, u_k, veh)

    x_k1 = x_k + (1/6)*dt*(w1 + 2*w2 + 2*w3 + w4)

    # adjust angles within bounds
    for d in eachindex(x_k1)
        if sg.angle_wrap_array[d] == true
            x_k1[d] = x_k1[d] % (2*pi)
            x_k1[d] > pi ? x_k1[d] -= 2*pi : x_k1[d] -= 0
        end
    end

    return x_k1
end

# workspace checker
function in_workspace(x, env, veh)
    veh_body = state_to_body(x, veh)

    if issubset(veh_body, env.workspace)
        return true
    end
        
    return false
end

# obstacle set checker
function in_obstacle_set(x, env, veh)
    veh_body = state_to_body(x, veh)
    
    for obstacle in env.obstacle_list
        if isdisjoint(veh_body, obstacle) == false
            return true
        end
    end

    return false
end

# target set checker
function in_target_set(x, env::Environment, veh::Vehicle)
    veh_body = state_to_body(x, veh)

    if issubset(veh_body, env.goal)
        return true
    end
        
    return false
end

# vehicle body transformation function
function state_to_body(x, veh)
    # rotate body about origin by theta
    rot_matrix = [cos(x[3]) -sin(x[3]); sin(x[3]) cos(x[3])]
    body = linear_map(rot_matrix, veh.origin_body)

    # translate body from origin by [x, y]
    trans_vec = x[1:2]
    LazySets.translate!(body, trans_vec)

    return body
end

function multi2single_ind(ind_m, sg)
    ind_s = 1
    for d in eachindex(ind_m)
        ind_s += (ind_m[d]-1)*prod(sg.state_grid.cut_counts[1:(d-1)])
    end

    return ind_s
end

function find_idx(val, array)
    idx = searchsortedfirst(array, val) - 1

    return idx
end

function plot_polygon(my_plot, P, lw, lc, ll)
    P_x_pts = [P[:,1]; P[1,1]]
    P_y_pts = [P[:,2]; P[1,2]]

    plot!(my_plot, P_x_pts, P_y_pts, linewidth=lw, linecolor=lc, label=ll)
end