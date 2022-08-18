
# obstacle set checker
function in_obstacle_set(y, O_set)
    for Oi_set in O_set
        rows = [collect(1:size(Oi_set, 1)); 1]

        ineqs = zeros(Bool, size(Oi_set, 1))
        for i in 1:size(Oi_set, 1)
            i1 = rows[i]
            i2 = rows[i+1]

            x1 = Oi_set[i1,1]
            y1 = Oi_set[i1,2]
            x2 = Oi_set[i2,1]
            y2 = Oi_set[i2,2]

            val = (y1 - y2)*y[1] + (x2 - x1)*y[2] + x1*y2 - x2*y1

            if val >= 0
                ineqs[i] = 1
            end
        end

        if all(ineqs) == true
            return true
        end
    end

    return false
end

# define obstacles
O1_set = [[-1.0 0.0];
        [1.0    0.0];
        [1.0    1.0];
        [-1.0   1.0]]

O2_set = [[1.5  -4.0];
        [2.0    -4.0];
        [2.0    -2.0];
        [1.5    -2.0]]

O_set = [O1_set, O2_set]

y = [1.5, -1, pi/2]

@show in_obstacle_set(y, O_set);

# ISSUE:
#   - returns everything as true