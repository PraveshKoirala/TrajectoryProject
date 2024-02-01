using Plots

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

function get_grid(knots, points=100)
    K = size(knots)[1]
    bins = get_bins(K)
    t = 0:(1/points):1
    return hcat([(1 .- t).^(K-n) .* t.^(n) .* bins[n+1] for n in 0:K]...)
end

function get_coords(knots, grid)
    println(size(grid))
    println(size(knots))
    return grid * vcat(knots, knots[1, :]')
end

function plot_splines(flags, knots, coords, xmin=-3, xmax=3, ymin=-3, ymax=3)
    p = plot(flags[:, 1], flags[:, 2], seriestype=:scatter, markershape=:utriangle, markercolor=:red, markersize=8,
        xlims=[xmin, xmax], ylims=[ymin,ymax])
    plot!(p, knots[:, 1], knots[:, 2], seriestype=:scatter, markercolor=:blue)
    plot!(p, coords[:, 1], coords[:, 2], ls=:dot)
    display(p)
end

function get_length(coords)
    # return length of bezier curve
    d = 0
    for i in 2:size(coords)[1]
        d += ((coords[i, 1] - coords[i-1, 1])^2 + (coords[i, 2]-coords[i-1,2])^2)^0.5
    end
    return d
end

flags = [-1 -1;
         -1 1;
         3 2;
         3 3]
knots1 = [
    -1 -1;
    2 1;
    3 3;
    3 4;
]

knots2 = [
    -2 -1;
    2 2;
    3 4;
    -3 4;
]

grid1 = get_grid(knots1)
grid2 = get_grid(knots2)
t1 = get_coords(knots1, grid1)
t2 = get_coords(knots2, grid2)
l1 = get_length(t1)
l2 = get_length(t2)

plot_splines(flags, knots1, get_coords(knots1, grid1))


# Two cars move in same speed, so we work with the one who has the shortest trajectory length
# and extrapolate for the longer one
short, shortlength, long, longlength = l1 < l2 ? t1, l1, t2, l2 : t2, l2, t1, l1

function collide(l1s, l1e, l2s, l2e)
    a1 = l1e[2] - l1s[2]
    b1 = l1s[1] - l1e[1]
    c1 = a1 * l1s[1] + b1 * l1[2]
    a2 = l2e[2] - l2s[2]
    b2 = l2s[1] - l2e[1]
    c2 = a2 * l2s[1] + b2 * l2s[2]
    Δ = a1 * b2 - a2 * b1
    intersect = [(b2 * c1 - b1 * c2) / Δ, (a1 * c2 - a2 * c1) / Δ]
    if min(l1s[1], l1e[1]) <= intersect[1] <= max(l1s[1], l1e[1]) &&
        min(l1s[2], l1e[2]) <= intersect[2] <= max(l1s[2], l1e[2]) &&
        min(l2s[1], l2e[1]) <= intersect[1] <= max(l2s[1], l2e[1]) &&
        min(l2s[2], l2e[2]) <= intersect[2] <= max(l2s[2], l2e[2]) &&
        return true, intersect  # intersects
    end
    return false, intersect # does not intersect
end

function score(flags, short, shortlength, long, longlength, fr=0.5)
    # short trajectory vs long trajectory
    T = size(short)[1]
    NF = size(flags)[1]
    ratio = shortlength / longlength
    
    captured = Dict([(i, 0) for i in 1:NF])
    distanceS = [Inf for _ in 1:NF]     # distance with all flags
    distanceL = [Inf for _ in 1:NF]     # distance with all flags

    ts = tl = 0
    last_long = long[0, :]
    last_short = short[0, :]
    while true
        ts += 1
        tl += ratio
        tl_a, tl_b = floor(tl), ceil(tl)
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
            if collides
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
            if captured[f] != 0 continue end
            ds = dl = false
            if short_coord != false
                ds = (flags[f, 1] - short_coord[1])^2 + (flags[f, 2] - short_coord[2])^2
                distanceS[f] = min(ds, distanceS[f])
            end
            dl = (flags[f, 1] - long_coord[1])^2 + (flags[f, 2] - long_coord[2])^2
            distanceL[f] = min(dl, distanceL[f])
            
            # check capture
            if ds != false && ds <= fr && dl <= fr
                # simultaneous capture!
                captured[f] = ds < dl ? -1 : 1
            elseif ds !=false && ds <= fr
                captured[f] = -1
            elseif dl <= fr
                captured[f] = 1
            end
        end

        # adjust ratio if short trajectory has ended
        if ts == T
            tl = floor(tl)
            ratio = 1
        end
        if tl >= T break end
    end

    #calculate score
    short_score = 0
    long_score = 0
    for f in 1:NF
        if captured[f] == -1
            short_score += flag_score
        elseif captured[f] == 1
            long_score += flag_score
        else
            short_score -= distance_factor * distanceS[f]
            long_score -= distance_factor * distanceL[f]
        end
    end
    return short_score, long_score
end
