# target set checker
function in_target_set(y, T_set)
    rows = [collect(1:size(T_set, 1)); 1]

    for i in 1:size(T_set, 1)
        i1 = rows[i]
        i2 = rows[i+1]

        x1 = T_set[i1,1]
        y1 = T_set[i1,2]
        x2 = T_set[i2,1]
        y2 = T_set[i2,2]

        val = (y1 - y2)*y[1] + (x2 - x1)*y[2] + x1*y2 - x2*y1

        println(val)

        if val < 0
            return false
        end
    end

    return true
end


# define target set
T_set = [[1.0   4.0];
        [2.0    4.0];
        [2.0    5.0];
        [1.0    5.0]]

y = [1.1, 5.1, pi/2]

@show in_target_set(y, T_set)