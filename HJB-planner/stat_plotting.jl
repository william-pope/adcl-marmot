using Plots

h_xy_arr = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.125, 0.1]

runtime1 = [2.968, 4.449, 4.875, 6.160, 7.448, 10.902, 16.676, 29.657, 43.459, 85.106, 166.976, 272.721, 604.331]
runtime2 = [0.0698, 0.1164, 0.1445, 0.1906, 0.2555, 0.4017, 0.6032, 1.175, 1.930, 3.451, 7.204, 12.347, 21.686]

p1 = plot(h_xy_arr, runtime2, xflip=true,
    xaxis=:log, yaxis=:log, size=(800,600),
    xticks=((0.1:0.1:1.0),["0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"]),
    # yticks=([2.968,43.459, 85.106, 166.976, 272.721, 604.331],string.([2.968,43.459, 85.106, 166.976, 272.721, 604.331])),
    yticks=([1.930, 3.451, 7.204, 12.347, 21.686],string.([1.930, 3.451, 7.204, 12.347, 21.686])),
    xlabel="XY Grid Spacing [m] (θ Spacing=10°)", ylabel="Time [s]", label="",
    title="HJB Calculation Time",
    markershape=:circle, markersize=7,
    linealpha=0)

plot!(p1, [0.25, 0.25], [minimum(runtime2), maximum(runtime2)],
    label="Approx Maximum Spacing", legend=:right,
    linecolor=:green, linestyle=:dash, linewidth=2)

# need minor ticks, would like axes in normal notation