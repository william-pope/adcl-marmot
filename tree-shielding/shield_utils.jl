# shield_utils.jl

function propagate_human(x_ih_k, ig, Dt, v_human, goal_positions)
    # break out current state
    xp_ih_k = x_ih_k[1]
    yp_ih_k = x_ih_k[2]

    # pull out chosen goal location
    xpg = goal_positions[ig].x
    ypg = goal_positions[ig].y

    # calculate derivative at current state
    C_x = ((xpg-xp_ih_k)^2 + (ypg-yp_ih_k)^2)^(-1/2)

    xp_ih_dot_k = v_human * C_x * (xpg-xp_ih_k)
    yp_ih_dot_k = v_human * C_x * (ypg-yp_ih_k)

    # calculate next state
    xp_ih_k1 = xp_ih_k + (xp_ih_dot_k * Dt)
    yp_ih_k1 = yp_ih_k + (yp_ih_dot_k * Dt)

    # reassemble state vector
    x_ih_k1 = [xp_ih_k1, yp_ih_k1]

    return x_ih_k1
end
