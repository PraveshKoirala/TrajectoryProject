using Serialization
using Plots
include("spline_methods.jl")


function plot_splines(flags, knots1, knots2, coords1, coords2, xmin=-10, xmax=10, ymin=-10, ymax=10)
    p = plot(flags[:, 1], flags[:, 2], seriestype=:scatter, markershape=:utriangle, markercolor=:red, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
    plot!(p, knots1[:, 1], knots1[:, 2], seriestype=:scatter, markercolor=:blue)
    plot!(p, knots2[:, 1], knots2[:, 2], seriestype=:scatter, markercolor=:red)

    c1 = cgrad(:blues, length(coords1), rev=true)
    c2 = cgrad(:heat, length(coords2), rev=true)
    for i in 1:size(coords1)[1]
        plot!(p, [coords1[i, 1]], [coords1[i, 2]], seriestype=:scatter, color=c1[i], legend=false, markersize=2, markerstrokewidth = 0)
        plot!(p, [coords2[i, 1]], [coords2[i, 2]], seriestype=:scatter, color=c2[i], legend=false, markersize=2, markerstrokewidth = 0)
    end
    
end

function plot_game(flags, coords1, coords2, xmin=-10, xmax=10, ymin=-10, ymax=10)
    d1, d2 = get_length(coords1), get_length(coords2)
    short, long = d1 < d2 ? (coords1, coords2) : (coords2, coords1)
    id1, fid1, id2, fid2 = d1 < d2 ? (:short_coord, -1, :long_coord, 1) : (:long_coord, 1, :short_coord, -1)
    history = score(flags, copy(short), get_length(short), copy(long), get_length(long), false, true)
    
    @gif for h in history
        
        p = plot(flags[:, 1], flags[:, 2], seriestype=:scatter, markershape=:utriangle, markercolor=:white, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
        for f in 1:size(flags)[1]
            fcolor = :white
            if h.captured[f][1] == fid1 
                fcolor = :blue 
            elseif h.captured[f][1] == fid2
                fcolor = :red 
            end
            p = plot!(p, [flags[f, 1]], [flags[f, 2]], seriestype=:scatter, markershape=:utriangle, markercolor=fcolor, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
        end
            
        plot!(p, coords1[:, 1], coords1[:, 2], ls=:dot, markercolor = :blue)
        plot!(p, coords2[:, 1], coords2[:, 2], ls=:dot, markercolor = :red)

        plot!(p, [h[id1][1]], [h[id1][2]], seriestype=:scatter, color=:blue, legend=false, alpha=0.5, markersize=12)
        plot!(p, [h[id2][1]], [h[id2][2]], seriestype=:scatter, color=:red, legend=false, alpha=0.5, markersize=12)
        plot!(p, [h[id1][1]], [h[id1][2]], seriestype=:scatter, color=:blue, legend=false)
        plot!(p, [h[id2][1]], [h[id2][2]], seriestype=:scatter, color=:red, legend=false)
        
    end
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

K = floor(Int64, length(KNOTS[end])/4)
knots1 = reshape(KNOTS[end][1:2*K], (2, K))'
knots2 = reshape(KNOTS[end][2*K+1:end], (2, K))'
grid1 = get_grid(knots1)
grid2 = get_grid(knots2)
t1 = get_coords(knots1, grid1)
t2 = get_coords(knots2, grid2)
plot_game(flags, t1, t2)