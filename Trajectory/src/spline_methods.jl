include("constants.jl")

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

function get_grid(knots, points=50)
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



function get_length(coords)
    # return length of bezier curve
    d = 0
    for i in 2:size(coords)[1]
        d += ((coords[i, 1] - coords[i-1, 1])^2 + (coords[i, 2]-coords[i-1,2])^2)^0.5
    end
    return d
end

function collide(l1e, l2e)
    if ((l1e[1] - l2e[1])^2 + (l1e[2] - l2e[2])^2)^0.5 < CR
        return true, l1e
    end
    return false, l1e
end


# Function to calculate distance between two points
function distance(p1, p2)
    return norm(p2 - p1)
end

function interpolate(a, b, d)
    D = distance(a, b)
    return a + (b - a) * d/D
end

# Function to interpolate and find points on the curve
function generate_uniform_points(curve :: Matrix{Float64}, velocity = 1)
    segment = 1
    NSegment = size(curve)[1]
    points = [curve[1, :]]
    new_point = Nothing
    while true
        remaining = velocity
        while true
            d = distance(curve[segment, :], curve[segment+1, :])
            if remaining < d
                new_point = interpolate(curve[segment, :], curve[segment+1, :], remaining)
                curve[segment, :] = new_point
                break
            end
            remaining -= d
            segment += 1
            if segment >= NSegment return points end
        end 
        push!(points, new_point)
    end
end

function score(flags, short, shortlength, long, longlength, throwcollision=true, simulate=false)
    # short trajectory vs long trajectory
    velocity = 0.2
    short = generate_uniform_points(short, velocity)
    long = generate_uniform_points(long, velocity)

    Ts = size(short)[1]
    Tl = size(long)[1]
    NF = size(flags)[1]
    
    captured = Dict([(i, (0, 0.0)) for i in 1:NF])
    distanceS = [(Inf, Inf) for _ in 1:NF]     # distance with all flags
    distanceL = [(Inf, Inf) for _ in 1:NF]     # distance with all flags

    # return if in simulation mode
    history = []
    ts = tl = 0

    while true
        ts += 1
        tl += 1
        short_coord = long_coord = Nothing
        if ts <= Ts
            short_coord = short[ts]
        else
            short_coord = false     # shorty completed trajectory
        end
        if tl <= Tl
            long_coord = long[tl]
        else
            long_coord = false     # shorty completed trajectory
        end
        if short_coord == false && long_coord == false break end
        if short_coord != false && long_coord != false
            # check collision of trajectories
            collides, collision_point = collide(short_coord, long_coord)
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

        # calculate score for flags
        for f in 1:NF
            if captured[f][1] != 0 continue end
            ds = dl = false
            if short_coord != false
                ds = ((flags[f, 1] - short_coord[1])^2 + (flags[f, 2] - short_coord[2])^2)^0.5
                distanceS[f] = min((ds, ts), distanceS[f])
            end
            if long_coord != false
                dl = ((flags[f, 1] - long_coord[1])^2 + (flags[f, 2] - long_coord[2])^2)^0.5
                distanceL[f] = min((dl, tl), distanceL[f])
            end
                        
            # check capture
            if ds != false && ds <= FR && dl !=false && dl <= FR
                # simultaneous capture!
                captured[f] = ds < dl ? (-1, ts) : (1, tl)
            elseif ds !=false && ds <= FR
                captured[f] = (-1, ts)
            elseif ds !=false && dl <= FR
                captured[f] = (1, tl)
            end
        end
        if simulate push!(history, (long_coord = copy(long_coord), short_coord=copy(short_coord == false ? short[Ts] : short_coord), 
                        captured=copy(captured))) end
    end

    #calculate score
    short_score = 0
    long_score = 0
    distance_factor = 0.1
    flag_score = 20
    for f in 1:NF
        if captured[f][1] == -1 # captured by short
            short_score += flag_score
            long_score -= distance_factor*distanceL[f][1]^0.5
        elseif captured[f][1] == 1  #captured by long
            long_score += flag_score
            short_score -= distance_factor*distanceS[f][1]^0.5
        else
            short_score -= distanceS[f][1]^0.5
            long_score -= distanceL[f][1]^0.5
        end
        
    end
    # println(captured)
    if simulate return history end
    return short_score - 0.1*shortlength, long_score -0.1*longlength
end