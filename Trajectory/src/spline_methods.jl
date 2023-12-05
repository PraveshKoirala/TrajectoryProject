using Plots
using Images

function get_bins(K=4)
    # For the given number of knots, find the binomial bezier coefficients. Hard-coded for now
    if K == 3
        return [1,3,3,1]
    elseif K == 4
        return [1, 4, 6, 4, 1]
    elseif K == 5
        return [1, 5, 10, 10, 5, 1]
    elseif K == 6
        return [1, 6, 15, 20, 15, 6, 1]
    end
    throw("Hardcoded coefficients only for 3,4,5 flags")
end

function get_grid(knots, points=20)
    K = size(knots)[1]
    bins = get_bins(K)
    t = 0:(1/points):1
    return hcat([(1 .- t).^(K-n) .* t.^(n) .* bins[n+1] for n in 0:K]...)
end

function get_coords(knots, grid)
    # println(size(grid))
    # println(size(knots))
    return grid * vcat(knots, knots[1, :]')
end

function plot_splines(flags, knots1, knots2, coords1, coords2, xmin=-10, xmax=10, ymin=-10, ymax=10)
    flagred = load("src/assets/flag.png")
    p = plot(flags[:, 1], flags[:, 2], seriestype=:scatter, markershape=:utriangle, markercolor=:red, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
    plot!(p, knots1[:, 1], knots1[:, 2], seriestype=:scatter, markercolor=:blue)
    plot!(p, knots2[:, 1], knots2[:, 2], seriestype=:scatter, markercolor=:red)
    plot!(p, coords1[:, 1], coords1[:, 2], ls=:dot)
    plot!(p, coords2[:, 1], coords2[:, 2], ls=:dot)
end

function get_length(coords)
    # return length of bezier curve
    d = 0
    for i in 2:size(coords)[1]
        d += ((coords[i, 1] - coords[i-1, 1])^2 + (coords[i, 2]-coords[i-1,2])^2)^0.5
    end
    return d
end

function collide(l1s, l1e, l2s, l2e)
    if ((l1e[1] - l2e[1])^2 + (l1e[2] - l2e[2])^2)^0.5 < 0.5
        return true, l1e
    end
    return false, l1e

    # a1 = l1e[2] - l1s[2]
    # b1 = l1s[1] - l1e[1]
    # c1 = a1 * l1s[1] + b1 * l1s[2]
    # a2 = l2e[2] - l2s[2]
    # b2 = l2s[1] - l2e[1]
    # c2 = a2 * l2s[1] + b2 * l2s[2]
    # Δ = a1 * b2 - a2 * b1
    # intersect = [(b2 * c1 - b1 * c2) / Δ, (a1 * c2 - a2 * c1) / Δ]
    # if min(l1s[1], l1e[1]) <= intersect[1] <= max(l1s[1], l1e[1]) &&
    #     min(l1s[2], l1e[2]) <= intersect[2] <= max(l1s[2], l1e[2]) &&
    #     min(l2s[1], l2e[1]) <= intersect[1] <= max(l2s[1], l2e[1]) &&
    #     min(l2s[2], l2e[2]) <= intersect[2] <= max(l2s[2], l2e[2]) &&
    #     return true, intersect  # intersects
    # end
    # return false, intersect # does not intersect
end

function score(flags, short, shortlength, long, longlength, 
    throwcollision=true, fr=0.5)
    # short trajectory vs long trajectory
    T = size(short)[1]
    NF = size(flags)[1]
    ratio = shortlength / longlength
    
    captured = Dict([(i, (0, 0.0)) for i in 1:NF])
    distanceS = [(Inf, Inf) for _ in 1:NF]     # distance with all flags
    distanceL = [(Inf, Inf) for _ in 1:NF]     # distance with all flags

    ts = tl = 1
    last_long = long[1, :]
    last_short = short[1, :]
    while true
        ts += 1
        tl += ratio
        tl_a, tl_b = floor(Int64, tl), ceil(Int64, tl)
        long_a = long[tl_a, :]
        long_b = long[tl_b, :]
        delta = long_b - long_a
        long_coord = long_a + (tl - tl_a) * delta
        if ts <= T
            short_coord = short[ts, :]
        else
            short_coord = false     # shorty completed trajectory
        end

        if short_coord != false
            # check collision of trajectories
            collides, collision_point = collide(last_short, short_coord, last_long, long_coord)
            if collides && throwcollision
                ds_c = (collision_point[1] - last_short[1])^2 + (collision_point[2] - last_short[2])^2
                dl_c = (collision_point[1] - last_long[1])^2 + (collision_point[2] - last_long[2])^2
                if ds_c < dl_c
                    # shorty reaches first
                    throw("Short collision")
                else
                    # longy reaches first
                    throw("Long collision")
                end
            end
        end
        last_short = short_coord
        last_long = long_coord

        # calculate score for flags
        for f in 1:NF
            if captured[f][1] != 0 continue end
            ds = dl = false
            if short_coord != false
                ds = ((flags[f, 1] - short_coord[1])^2 + (flags[f, 2] - short_coord[2])^2)^0.5
                distanceS[f] = min((ds, ts), distanceS[f])
            end
            dl = ((flags[f, 1] - long_coord[1])^2 + (flags[f, 2] - long_coord[2])^2)^0.5
            distanceL[f] = min((dl, tl), distanceL[f])
            
            # check capture
            if ds != false && ds <= fr && dl <= fr
                # simultaneous capture!
                captured[f] = ds < dl ? (-1, ts) : (1, tl)
            elseif ds !=false && ds <= fr
                captured[f] = (-1, ts)
            elseif dl <= fr
                captured[f] = (1, tl)
            end
        end

        # adjust ratio if short trajectory has ended
        if ts == T
            tl = floor(Int64, tl)
            ratio = 1
        end
        if tl >= T break end
    end

    #calculate score
    short_score = 0
    long_score = 0
    distance_factor = 5
    flag_score = 30
    for f in 1:NF
        if captured[f][1] == -1 # captured by short
            short_score += flag_score
            long_score -= distanceL[f][1]^0.5
        elseif captured[f][1] == 1  #captured by long
            long_score += flag_score
            short_score -= distanceS[f][1]^0.5
        else
            short_score -= distance_factor*distanceS[f][1]^0.5
            long_score -= distance_factor*distanceL[f][1]^0.5
        end
        
    end
    return short_score, long_score
end
