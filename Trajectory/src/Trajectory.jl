# Uses solver to extract a feasible point

using Symbolics
using LinearAlgebra
using Random

Random.seed!(20)

if !(@isdefined POLICY)
    println("Please specify which policy to use.")
    println("Set POLICY=1 for linear and POLICY=2 for non-linear.")
    throw(ArgumentError("Policy not set."))
end

if POLICY == 1
    println("Working with Linear Policy. An initial feasible point will be obtained using the Solver.")
    include("trajectory_tri_level.jl")
else
    println("Working with NonLinear Policy. An initial feasible point is being used manually.")
    include("trajectory_tri_level_nonlin.jl")
end

# Gradient of the upper function w.r.t all decision variables.
JF_tau = Symbolics.gradient(F_tau(x_), x_)
JF_t = Symbolics.gradient(F_t(x_), x_)
JF_T = Symbolics.gradient(F_T(x_), x_)

JF = [JF_tau, JF_t, JF_T]

JG_t = Symbolics.jacobian(G_t(x_), x_)
JG_T = Symbolics.jacobian(G_T(x_), x_)

JG = [[], JG_t, JG_T]

function evaluate(expression, x)
    vars_mapping = Dict(zip(x_, x))
    return substitute(expression, vars_mapping)
end


function is_feasible(x, I=vec(1:G_length), ϵ=1e-6)
    # return true if a feasible point
    # By default, calculate for every constraint.
    
    GI = []
    for g in G[I]
        push!(GI, evaluate(g, x))
    end

    GI = [gi.val for gi in GI]  # convert to float. MUST BE A BETTER WAY!
    # println(GI)

    if any(GI .< -ϵ)
        # if point is infeasible, return false and the calculated constraint violations as well
        return false, GI
    end
    return true, GI
end

function get_upper_moves(N, x, player_level, maintain_feasibility=true)
    # N: Number of random samples to generate
    # x: The variable vector
    
    I = decision_variables[player_level]
    Jac = JF[player_level]

    # must satisfy constraints of *all* bottom players
    cset = constraints_map[player_level]

    # I: The index set that this player controls
    # Jac: Jacobian of this function

    alpha = 2
    beta = 0.1
    
    # one of the directions returned is player 1's gradient
    x_g = -evaluate(Jac, x) .* beta
    x_g = [t.val for t in x_g]   # convert to float. must be a better way.

    d = length(x)
    num_directions = 0
    
    Ī = setdiff(Set(1:d), I)    # the variables that this player does NOT control
    Ī = [_ for _ in Ī]  # there should be an easy conversion to be honest
    
    # list of candidate directions
    x_directions = []
    # direction of the gradient    
    new_direction = []
    infeasible_attempts = 0
    while num_directions < N
        new_direction = isempty(new_direction) ? x_g : (rand(d) .- 0.5 ) .* alpha
        new_direction[Ī] .= 0
        x_new = x + new_direction
        if maintain_feasibility
            feasible, slack = is_feasible(x_new, cset)
            if !feasible 
                infeasible_attempts += 1
                # Tried more than 2N times but couldn't get a good point,
                if infeasible_attempts > 2*N
                    return isempty(x_directions) ? [] : x .+ x_directions
                end
                continue
            end
        end
        x_directions = isempty(x_directions) ? x_new : hcat(x_directions, new_direction)
        num_directions += 1
    end
    return x.+x_directions
end


function get_lower_moves(x, I, Jac, player_level=3)
    # x: Variable vector
    # I: index set of variables that this player controls
    # Jacobian of this player

    # move towards the -ve gradient while maintaining lower feasibility
    alpha = 1
    # also return stuffs like norm zero etc.

    d=length(x)
    Ī = setdiff(Set(1:d), I)    # the variables that this player does NOT control
    Ī = [_ for _ in Ī]  # there should be an easy conversion to be honest

    grad = function (x)
        sd = -evaluate(Jac, x)
        # sd is a Num array so convert to float and select only controllable gradients
        sd = [s.val for s in sd]
        sd[Ī] .= 0
        return sd
    end
    
    c_set = constraints_map[player_level]
    
    # Make this point feasible
    # This leverages the affine property of constraints 
    # For general, this won't work
    for _ = 1:100
        # initial feasibility
        feasible, slack = is_feasible(x, c_set)
        if feasible break end
        most_violated_constraint = argmin(slack)
        # direction towards the constraint
        jd = JG[player_level][most_violated_constraint, I]
        jd = evaluate(jd, x)
        x[I] += [j.val for j in jd]
    end
    
    for _ = 1:20
        sd = grad(x)
        N = abs(norm(sd))
        if (N < 1e-4)
            break
        end

        for __ = 1:10
            x_n = x + alpha*sd
            feasible, _ = is_feasible(x_n, c_set, 1e-6)    #check feasibility of just the bottom player
            # Norm of gradient at new point should consistently decrease
            new_N = abs(norm(grad(x_n)))
            if feasible && (new_N <= N)
                x = x_n
                @goto OUTER_LOOP
            end
            # Norm didn't decrease (or entered infeasible region), so adjust step size
            alpha /= 10
        end
        # no solution found, return
        
        break
        @label OUTER_LOOP
        if alpha < 1e-4
            # Step size too little, no substantial improvement possible in future
            break
        end
    end
    # println("Leaving with player2 norm=", N)
    return x
end

# Uses solver to find feasible point inside the constraint boundary
println("Starting with a feasible point found for the problem: ", x_s)
# Now that we have a feasible point, we can execute the montecarlo method.

function get_next_step(n, x, player_id)
    I = decision_variables[player_id]
    Jf = JF[player_id]
    f = F[player_id]

    if player_id == L
        # This is the bottom level player
        x = get_lower_moves(x, I, Jf, player_id)
        return x
    end

    # This is AD-HOC. Need to figure out a better way of selecting optimal points
    f_val_2 = f(x)

    c_set = constraints_map[player_id]
    # does this point satisfy constraint? (It must when starting out)
    feasible, _ = is_feasible(x, c_set)

    X_t = get_upper_moves(n, x, player_id)
    X_t = isempty(X_t) ? x[:,:] : hcat(x, X_t)
    num_points = size(X_t)[2]
    candidates = []
    for m = 1:num_points
        x_t = X_t[:, m]
        x_t = get_next_step(n, x_t, player_id+1)
        if x_t == Nothing
            # couldn't find a feasible point
            continue
        end
        f_xt = f(x_t)
        
        feasible, _ = is_feasible(x_t, c_set)
        # This is ad-hoc. Need a better way.
        if feasible && (player_id == 2 ? f_xt <= f_val_2 : true)
            push!(candidates, (f_xt, x_t))
        end
    end
    if length(candidates) == 0
        # couldn't find any feasible direction
        return Nothing
    end
    min_fx, min_x = minimum(candidates)
    # Best possible result you could get from here
    return min_x
end

# number of gradient to sample for top player
if !(@isdefined N)
    N = 10
end
println("Running with samples per player N = ", N, ". To change it, just set the value in the terminal.")
println()

didnt_update_since = 0
if !(@isdefined MAX_ITER)
    MAX_ITER = 20
end

println("Running with MAX_ITER = ", MAX_ITER, ". To change it, just set the value in the terminal.")
# path
px = [x_s[1]]
py = [x_s[2]]
for i = 1:MAX_ITER
    global x_s, didnt_update_since, px, py
    println("Iteration: ", i)
    if didnt_update_since > 10
        # Stagnated
        println("Stagnated, breaking now...")
        break
    end

    last = x_s
    x_s = get_next_step(N, x_s, 1)

    if x_s == Nothing
        println("No new points were found. Keeping this one.")
        x_s = last
    end

    if all(x_s ≈ last)
        println("Continuing with the same point.")
        didnt_update_since += 1
    else
        println("Better point found in the neighborhood.")
        println(x_s)
        push!(px, x_s[1])
        push!(py, x_s[2])
        println("Old objective: ", F_tau(last), ". New objective: ", F_tau(x_s))
        println()
        println()
        didnt_update_since = 0
    end
end

println("Concluded with the following statistics:")
println("Top objective= ", evaluate(F[1](x_s), x_s))
println("x= ", x_s)
println("Feasible?= ", is_feasible(x_s)[1])

X = x_s[xdim]
function circleShape(h, k, r)
    θ = LinRange(0, 2*π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

using Plots

τx = []
τy = []
tx, ty = X[1], X[2]
for n in 1:N_τ
    global tx, ty, τx, τy
    τx = push!(τx, tx)
    τy = push!(τy, ty)
    tx, ty = PI([tx, ty])
end

# feasible region
plot(circleShape(r, r, r), seriestype=[:shape,], c=:blue, linecolor= :black, fillalpha=0.2, aspect_ratio=1)
#obstacle
plot!(circleShape(O[1], O[2], or), seriestype=[:shape,], c=:red, fillalpha=0.2, aspect_ratio=1)
#value of x
plot!([X[1]], [X[2]], seriestype=:scatter, c=:blue)

#Draw x1=D plane
plot!([D for i in 0:20], 0:20)

#Draw taus
plot!(LinRange(X[1], D, 20), [X[2] for i in 1:20], c=:blue, ls=:dash)
plot!(τx, τy, c=:purple, seriestype=:scatter)
plot!(px, py, c=:red, seriestype=:path)

