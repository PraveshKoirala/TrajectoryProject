using Serialization
using Plots
include("spline_methods.jl")

flags = [
0 2;
2 2;]

function plot_splines(flags, knots1, knots2, coords1, coords2, xmin=-10, xmax=10, ymin=-10, ymax=10)
    flagred = load("src/assets/flag.png")
    p = plot(flags[:, 1], flags[:, 2], seriestype=:scatter, markershape=:utriangle, markercolor=:red, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
    plot!(p, knots1[:, 1], knots1[:, 2], seriestype=:scatter, markercolor=:blue)
    plot!(p, knots2[:, 1], knots2[:, 2], seriestype=:scatter, markercolor=:red)
    plot!(p, coords1[:, 1], coords1[:, 2], ls=:dot)
    plot!(p, coords2[:, 1], coords2[:, 2], ls=:dot)
end

function plot_game(flags, coords1, coords2, k=1, xmin=-10, xmax=10, ymin=-10, ymax=10)
    flagred = load("src/assets/flag.png")
    p = plot(flags[:, 1], flags[:, 2], seriestype=:scatter, markershape=:utriangle, markercolor=:red, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
    k1 = k2 = k
    if size(coords1)
    plot!(p, coords1[k, 1], coords1[k, 2], seriestype=:scatter, markercolor=:blue)
    plot!(p, coords2[k, 1], coords2[k, 2], seriestype=:scatter, markercolor=:red)
    plot!(p, coords1[:, 1], coords1[:, 2], ls=:dot)
    plot!(p, coords2[:, 1], coords2[:, 2], ls=:dot)
end

KNOTS = deserialize("trajectories.dat")
@gif for k = 1:size(KNOTS)[1]
    knots = KNOTS[k]
    K = floor(Int64, length(knots)/4)
    knots1 = reshape(knots[1:2*K], (2, K))'
    knots2 = reshape(knots[2*K+1:end], (2, K))'

    grid1 = get_grid(knots1)
    grid2 = get_grid(knots2)
    t1 = get_coords(knots1, grid1)
    t2 = get_coords(knots2, grid2)

    plot_splines(flags, knots1, knots2, t1, t2)
end

@gif for k=1: