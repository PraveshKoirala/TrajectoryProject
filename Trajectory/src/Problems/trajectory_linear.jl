module LinearTrajectory
    include("Problems.jl")
    using .Problems
    using Plots
    using JuMP
    using Ipopt
    
    # Setup the environment
    r = 5        #radius of the feasibility region in positive quadrant
    R = [r, r]  # center of the feasibility circle

    O = [12, 9]  # obstacle center
    or = 2      # radius of the obstacle

    N_τ = 20     #number of trajectory points is n+1
    D = 20      #destination plane, any value larger than the obstacle should be fine

    
    function LinearPI(xnm1)
        δ = 1
        # Linear policy, just move δ step forward in x direction
        return xnm1[1] + δ, xnm1[2]
    end

    PI = LinearPI

    # Tri-level problem
    L = 3

    dim = 4;
    xdim = 1:2;
    tdim = 3:3;
    Tdim = 4:4;
    
    # Slack from the obstacle. >=0 if trajectory is not touching the obstacle.
    function g_tau(tau)
        return (tau[1]-O[1])^2 + (tau[2]-O[2])^2 - or^2
    end

    function chi(x0)
        # check if the initial point in feasible region
        return r^2-(x0[1]-R[1])^2-(x0[2]-R[2])^2
    end

    function F_T(x)
        return -x[Tdim][1]
    end

    function G_T(x)
        T = x[Tdim][1]
        constraints = [T, g_tau(x) - T];
        x_next = x;
        # Generate next trajectory points using the policy
        for _ = 2:N_τ
            x_next = PI(x_next)
            push!(constraints, g_tau(x_next) - T)
        end
        return constraints
    end

    # objective for upper
    function F_t(x)
        t = x[tdim][1]
        return .5*t^2
    end

    function G_t(x)
        T = x[Tdim][1];
        t = x[tdim][1]
        return [
            t-T,
            chi(x),
            t
        ]
    end


    function F_tau(x)
        return - (D - x[1])
    end

    function G_tau(x)
        return [1]
    end
    
    # four variables
    TrajectoryProblem = MultiLevelProblem(4)
    
    # For the first level, the second objective is more important because of the shared degree of freedom
    TrajectoryProblem.addLevel!(Problem((F_t, F_tau), G_tau, xdim))
    TrajectoryProblem.addLevel!(Problem(F_t, G_t, tdim))
    TrajectoryProblem.addLevel!(Problem(F_T, G_T, Tdim))
    # TrajectoryProblem.visualize = visualize

    function TrajectoryProblem.visualize(x; kwargs...)
        X = x[1:2]
        function circleShape(h, k, r)
            θ = LinRange(0, 2*π, 500)
            h .+ r*sin.(θ), k .+ r*cos.(θ)
        end


        τx = []
        τy = []
        tx, ty = X[1], X[2]
        for n in 1:N_τ
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
        px, py = kwargs[:px], kwargs[:py]
        plot!(px, py, c=:red, seriestype=:path)
    end

    function find_feasible_point(x_s)
        # Use Ipopt to find an initial feasible point.
        # The point will be inside the constraint boundary but not exactly feasible in Multilevel optimization sense
        # After a point is obtained, we must let the lower level players optimize to truly obtain a feasible starting point.
        
        # To find a starting point, we linearly weigh the function objectives with the last level getting the most
        # weight. This is in hopes that we'll achieve something sensible or near the local optimum.
        r1 = 1;
        r2 = 10;
        r3 = 100;
    
        m = Model(Ipopt.Optimizer)
        @variable(m, x[i=1:dim])
    
        C = vcat(G_T(x), G_t(x))
        lc = length(C)
    
        @NLconstraint(m, [i=1:lc], C[i]>=0)
    
        obj = @expression(m, r1 * F_tau(x) + r2 * F_t(x) + r3 * F_T(x));
        @NLobjective(m, Min, obj)
        optimize!(m)
        println(termination_status(m))
        return value.(x)
    end

    TrajectoryProblem.x_s = [0.015401027606127604, 5.14272294205888, 11, 10.878760000000172];
    # TrajectoryProblem.x_s = [1.5702263167931156, 3.5173743496664365, 26.25196450419984, 26.24388944799616]
    TrajectoryProblem.x_s = find_feasible_point([0, 0])

    export TrajectoryProblem
end