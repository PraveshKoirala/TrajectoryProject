module Trajectory
    # include("solvefull.jl")
    # include("Problems/trajectory_linear.jl")
    # include("Problems/trajectory_nonlinear.jl")
    include("Problems/AnandalingamDecentralized.jl")
    
    # include("Problems/Sinha.jl")
    # include("Problems/Tilahun.jl")
    # include("Problems/NestedToll.jl")

    # using Symbolics
    using LinearAlgebra
    using Random
    using JuMP
    using Statistics
    import Symbolics
    using Ipopt
    using Base.Threads

    # p for problem
    # using .LinearTrajectory: TrajectoryProblem as p
    # using .NonLinearTrajectory: TrajectoryProblem as p
    using .ALEx: ALProblem as p
    # using .SinhaEx1: SinhaProblem as p
    # using .Tilahun: TilahunProblem as p
    # using .NToll: NTollProblem as p


    function evaluate(expression, x)
        vars_mapping = Dict(zip(p.X, x))
        return [v.val for v in Symbolics.substitute(expression, vars_mapping)]
    end


    function is_feasible(x, I=vec(1:p.levels()), ϵ=1e-6)
        # return true if a feasible point
        # By default, calculate for every player.
        for i in I
            if any(p[i].g(x) .< -ϵ) return false end
        end
        return true
    end

    function get_upper_moves(N, x, player_level, maintain_feasibility=false)
        # N: Number of random samples to generate
        # x: The variable vector
        
        # I: The index set that this player controls
        # Jac: Jacobian of this function
        I = p[player_level].I
        # Jac = p.Jf(player_level)       # jacobian of this player

        # must satisfy constraints of *all* bottom players
        # cset = player_level:p.levels()

        alpha = p.alpha
        
        # one of the directions returned is this player's gradient
        # x_g = -evaluate(Jac, x) .* beta
        d = length(x)
        num_directions = 0
        
        Ī = setdiff(Set(1:d), I)    # the variables that this player does NOT control
        Ī = [_ for _ in Ī]  # there should be an easy conversion to be honest
        
        # list of candidate directions
        x_directions = []
        # direction of the gradient    
        new_direction = []
        while num_directions < N
            # new_direction = isempty(new_direction) ? x_g : (rand(d) .- 0.5 ) .* alpha
            new_direction = (rand(d) .- 0.5 ) .* alpha
            new_direction[Ī] .= 0
            x_directions = isempty(x_directions) ? new_direction : hcat(x_directions, new_direction)
            num_directions += 1
        end
        return x.+x_directions
    end

    function get_lower_moves(x_val, I, player_level=3)
        # Eventually we want to solve the final level using a solver like PATH but for testing purposes,
        # let's encode the final complete solution as high alpha and large number of samples.
        return x_val
        # x: Variable vector
        # I: index set of variables that this player controls
        # Jacobian of this player (not needed anymore)
        # dim = length(p.X)
        
        # @Symbolics.variables x[1:dim]
        # # Get objective and constraint of the lower:
        # O = p[player_level].f(x)[end]
        # C = p[player_level].g(x)
        
        # global m = Nothing
        # global x = Nothing
        
        # m = Model(Ipopt.Optimizer)
        # # set_attribute(m, "max_iter", 5)
        # set_silent(m)

        # # todo: use the start to start closer
        # @variable(m, x[i=1:dim])

        # # player only controls their own variables
        # Ī = setdiff(1:dim, I)
        # @constraint(m, x[Ī] .== x_val[Ī])
        
        # # Add constraint
        # for i in eachindex(C)
        #     # println("@NLconstraint(m, $(repr(C[i])) >= 0)")
        #     eval(Meta.parse("@NLconstraint(m, $(repr(C[i])) >= 0)"))
        # end
        # # Add objective
        # eval(Meta.parse("@NLobjective(m, Min, $(repr(O)))"))

        # optimize!(m)
        # return value.(x)
    end

    function get_next_step(n, x, player_id)
        n = p[player_id].n
        I = p[player_id].I
        # Jf = p.Jf(player_id)
        f = p[player_id].f
        # c_set = player_id:p.levels()
        c_set = player_id:player_id

        if player_id == p.levels()
            # This is the bottom level player
            x = get_lower_moves(x, I, player_id)
            return is_feasible(x, c_set) ? x : Nothing
        end

        # does this point satisfy constraint? (It must when starting out for all lower)

        X_t = get_upper_moves(n, x, player_id)
        
        num_points = size(X_t)
        if length(num_points) == 1
            if num_points[0] == 0 return Nothing end    # couldn't sample anything
            num_points = 1
        else
            num_points = num_points[2]
        end

        X_t = hcat(x, X_t)  # also include the zero direction
        candidates = []
        for m = 1:num_points
            x_t = X_t[:, m]
            x_t = get_next_step(n, x_t, player_id+1)
            if x_t == Nothing
                # couldn't find a feasible point
                continue
            end
            f_xt = f(x_t)
            
            feasible = is_feasible(x_t, c_set)
            if feasible
                push!(candidates, (f_xt, x_t))
            end
        end
        if length(candidates) == 0
            # couldn't find any feasible direction
            return Nothing
        end
        f_values = [c[1] for c in candidates]   # tuple of all function f_values
        # for a single function in each level, just return the minimizer
        if length(f_values[1]) == 1 
            min_fx, min_x = minimum(candidates)
            return min_x
        end
        # to matrix
        f_values = hcat([collect(e) for e in f_values]...)'
        mins = minimum(f_values, dims=1)
        maxes = maximum(f_values, dims=1)
        # construct μ for each fᵢ that maps minimum fᵢ-> 1 and maximum fᵢ-> ϵ will represent willingness to move
        ϵ = 0.1     #can't let mu be 0
        mu(f_v) = max.(1 .+ (f_v .- mins)./(mins.-maxes), zero(f_v).+ϵ)
        d(f_v) = prod(mu(f_v), dims=2)
        favorable_idx = argmax(d(f_values))[1]
        return candidates[favorable_idx][2]
    end

    function approximate(P, K=10)
       K = min(K, length(P))
       P = P[end-K+1:end]   # last k samples
       f = map(P) do xk
            map(i -> p[i].f(xk)[1], 1:p.levels())
       end
       return P[argmin(f)]
    end

    function run_opt()
        x_s = p.x_s
        # approximate smoothing
        P = []
        # Uses solver to find feasible point inside the constraint boundary
        println("Starting with a feasible point found for the problem: ", x_s)
        # Now that we have a feasible point, we can execute the montecarlo method.

        # number of gradient to sample for top player
        if !(@isdefined N)
            N = 5
        end
        println("Running with samples per player N = ", N, ". To change it, just set the value in the terminal.")
        println()

        didnt_update_since = 0
        MAX_ITER = p.MAX_ITER

        println("Running with MAX_ITER = ", MAX_ITER, ". To change it, just set the value in the terminal.")
        # path
        px = [x_s[1]]
        py = [x_s[2]]
        for i = 1:MAX_ITER
            # global x_s, didnt_update_since, px, py
            println("Iteration: ", i)
            if didnt_update_since > 20
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
            push!(P, x_s)
            println(x_s)
            if all(x_s ≈ last)
                println("Continuing with the same point.")
                didnt_update_since += 1
            else
                println("Better point found in the neighborhood.")
                println(x_s)
                push!(px, x_s[1])
                push!(py, x_s[2])
                println("Old objective: ", p[1].f(last), ". New objective: ", p[1].f(x_s))
                println()
                println()
                didnt_update_since = 0
            end
            # Decrease alpha per iteration
            p.alpha /= p.cooldown
        end
        
        println("Concluded with the following statistics:")
        println("Top objective= ", p[1].f(x_s))
        println("x= ", x_s)
        println("Feasible?= ", is_feasible(x_s))

        println("Smoothed solution:")
        # Remove this line to remove the smoothing
        x_s = approximate(P)
        println("Top objective= ", p[1].f(x_s))
        println("x= ", x_s)
        println("Feasible?= ", is_feasible(x_s))

        p.visualize(x_s; px=px, py=py)
    end

    export run_opt
    export p
end
