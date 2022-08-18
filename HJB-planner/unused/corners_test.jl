
using Plots

struct Vehicle
    c_vf::Float64   # max forward speed [m/s]
    c_vb::Float64   # max backward speed [m/s]
    c_phi::Float64  # max steering angle [rad]
    wb::Float64     # wheelbase [m]
    l::Float64      # length [m]
    w::Float64      # width [m]
    b2a::Float64    # rear bumber to rear axle [m]
end

function pose_to_corners(y, veh::Vehicle)
    x_BR = y[1] - veh.b2a*cos(y[3]) + (veh.w/2)*sin(y[3])
    x_BL = y[1] - veh.b2a*cos(y[3]) - (veh.w/2)*sin(y[3])
    x_FR = y[1] + (veh.l - veh.b2a)*cos(y[3]) + (veh.w/2)*sin(y[3])
    x_FL = y[1] + (veh.l - veh.b2a)*cos(y[3]) - (veh.w/2)*sin(y[3])

    y_BR = y[2] - veh.b2a*sin(y[3]) - (veh.w/2)*cos(y[3])
    y_BL = y[2] - veh.b2a*sin(y[3]) + (veh.w/2)*cos(y[3])
    y_FR = y[2] + (veh.l - veh.b2a)*sin(y[3]) - (veh.w/2)*cos(y[3])
    y_FL = y[2] + (veh.l - veh.b2a)*sin(y[3]) + (veh.w/2)*cos(y[3])
    
    veh_corners = [[x_BL, y_BL],
                    [x_BR, y_BR],
                    [x_FR, y_FR],
                    [x_FL, y_FL]]

    return veh_corners
end

marmot = Vehicle(1.5, 0.75, 0.475, 0.324, 0.5207, 0.2762, 0.0889)   

y = [0, 0, -pi/6]

veh_corners = pose_to_corners(y, marmot)

display(veh_corners)

p_corners = plot(aspect_ratio=:equal)
plot!(p_corners, getindex.(veh_corners,1), getindex.(veh_corners,2), markershape=:circle, markersize=4)
plot!(p_corners, [y[1]], [y[2]], markershape=:circle, markersize=4)