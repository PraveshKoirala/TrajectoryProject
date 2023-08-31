module LinearTrajectory
    include("Problems.jl")
    using .Problems
    using Plots
    using JuMP
    using Ipopt
    
    # Setup the environment
    r = 5        #radius of the feasibility region in positive quadrant
    R = [r, r]  # center of the feasibility circle

    # O = [15, 13]  # Top objective= (-17.350494166717795,)
                    # x= [2.6495058332822037, 9.346829770275754, 10.135786580949945, 9.46849889575728]
                    
    # O = [12, 9]   # Top objective= (-19.592984358059056,)
                    # x= [0.4070156419409441, 6.955217184755138, 0.40478628047151544, 0.34679850181433547]
                    
    # O = [15, 5]   # Top objective= (-19.471897308511,)
                    # x= [0.5281026914889992, 3.056468958195349, 0.02871955068917198, -9.880900609996543e-9]
                    
    O = [12, -3]    # Top objective= (-16.715481849682547,)
                    # x= [3.2845181503174548, 0.46271391595835354, 9.71352702589218, 8.071338249138893]
                    
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
        return t
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
    TrajectoryProblem.addLevel!(Problem(F_tau, G_tau, xdim, 2))
    TrajectoryProblem.addLevel!(Problem(F_t, G_t, vcat(xdim, tdim), 10))
    TrajectoryProblem.addLevel!(Problem(F_T, G_T, Tdim, 1))
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
            if tx > D
                break
            end
        end

        # feasible region
        plot(circleShape(r, r, r), label="Feasible Region", seriestype=[:shape,], c=:blue, linecolor= :black, fillalpha=0.2, aspect_ratio=1)
        #obstacle
        plot!(circleShape(O[1], O[2], or), label="Obstacle", seriestype=[:shape,], c=:red, fillalpha=0.2, aspect_ratio=1)
        

        #Draw x1=D plane
        plot!([D for i in 0:20], 0:20, label="Finish Line")

        #Draw taus
        # plot!(LinRange(X[1], D, 20), [X[2] for i in 1:20], c=:blue, ls=:dash)
        plot!(τx, τy, c=:purple, seriestype=:scatter, label="Trajectory")
        px, py = kwargs[:px], kwargs[:py]
        plot!(px, py, c=:red, seriestype=:path, label="Path", alpha=0.2)
        #value of x
        plot!([X[1]], [X[2]], label="Initial Point", seriestype=:scatter, c=:blue)
    end

    function find_feasible_point(x_s)
        # Use Ipopt to find an initial feasible point.
        # The point will be inside the constraint boundary but not exactly feasible in Multilevel optimization sense
        # After a point is obtained, we must let the lower level players optimize to truly obtain a feasible starting point.
        
        # To find a starting point, we linearly weigh the function objectives with the last level getting the most
        # weight. This is in hopes that we'll achieve something sensible or near the local optimum.
        r1 = 1e5;
        r2 = 1e-5;
        r3 = 1;
    
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

    # TrajectoryProblem.x_s = [0.015401027606127604, 5.14272294205888, 11, 10.878760000000172];
    # TrajectoryProblem.x_s = [1.5702263167931156, 3.5173743496664365, 26.25196450419984, 26.24388944799616]
    TrajectoryProblem.x_s = find_feasible_point([0, 0])
    TrajectoryProblem.MAX_ITER = 150
    TrajectoryProblem.alpha = 3
    export TrajectoryProblem
end