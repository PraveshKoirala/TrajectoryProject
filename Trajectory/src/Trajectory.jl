module Trajectory
    # include("Problems/trajectory_linear.jl")
    # include("Problems/trajectory_nonlinear.jl")
    include("Problems/AiyoshiShimizu1984Ex2.jl")
    include("Problems/Bard1988Ex2.jl")

    # using Symbolics
    using LinearAlgebra
    using Random
    using JuMP
    import Symbolics
    using Ipopt
    # p for problem
    # using .Ridge: RidgeProblem as p
    # using .LinearTrajectory: TrajectoryProblem as p
    # using .NonLinearTrajectory: TrajectoryProblem as p
    # using .AS1984: AS1984Problem as p
    using .Bard1988Ex2: Bard1988Ex2Problem as p

    Random.seed!(20)

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
        Jac = p.Jf(player_level)       # jacobian of this player

        # must satisfy constraints of *all* bottom players
        cset = player_level:p.levels()

        alpha = 2
        beta = 0.1
        
        # one of the directions returned is this player's gradient
        x_g = -evaluate(Jac, x) .* beta

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
                feasible = is_feasible(x_new, cset)
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

    function get_lower_moves(x_val, I, Jac, player_level=3)
        # x: Variable vector
        # I: index set of variables that this player controls
        # Jacobian of this player (not needed anymore)
        dim = length(p.X)
        
        @Symbolics.variables x[1:dim]
        # Get objective and constraint of the lower:
        O = p[player_level].f(x)[end]
        C = p[player_level].g(x)
        
        global m = Nothing
        global x = Nothing
        
        m = Model(Ipopt.Optimizer)
        set_silent(m)

        # todo: use the start to start closer
        @variable(m, x[i=1:dim])

        # player only controls their own variables
        Ī = setdiff(1:dim, I)
        @constraint(m, x[Ī] .== x_val[Ī])
        
        # Add constraint
        for i in eachindex(C)
            eval(Meta.parse("@NLconstraint(m, $(repr(C[i])) >= 0)"))
        end
        # Add objective
        eval(Meta.parse("@NLobjective(m, Min, $(repr(O)))"))

        optimize!(m)
        return value.(x)
    end

    function get_next_step(n, x, player_id)
        I = p[player_id].I
        Jf = p.Jf(player_id)
        f = p[player_id].f

        if player_id == p.levels()
            # This is the bottom level player
            x = get_lower_moves(x, I, Jf, player_id)
            return x
        end

        c_set = player_id:p.levels()
        # does this point satisfy constraint? (It must when starting out for all lower)
        feasible = is_feasible(x, c_set)

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
            
            feasible = is_feasible(x_t, c_set)
            if feasible
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

    function run_opt()
        x_s = p.x_s
        # Uses solver to find feasible point inside the constraint boundary
        println("Starting with a feasible point found for the problem: ", x_s)
        # Now that we have a feasible point, we can execute the montecarlo method.

        # number of gradient to sample for top player
        if !(@isdefined N)
            N = 10
        end
        println("Running with samples per player N = ", N, ". To change it, just set the value in the terminal.")
        println()

        didnt_update_since = 0
        if !(@isdefined MAX_ITER)
            MAX_ITER = 100
        end

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
        end

        println("Concluded with the following statistics:")
        println("Top objective= ", p[1].f(x_s))
        println("x= ", x_s)
        println("Feasible?= ", is_feasible(x_s))

        p.visualize(x_s; px=px, py=py)
    end

    export run_opt
end
