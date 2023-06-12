module Ridge
    include("Problems.jl")
    using .Problems

    w = (1, 2)
    function F(X)
        return 0.5 * (X[1]-w[1])^2 + 0.5* (X[2]-w[2])^2
    end

    # Unconstrained
    function G(X)
        return []
    end
    
    I = [1, 2]      # upper controls both x1 and x2

    function f(X)
        return 0.5 * (X[1]-X[2])^2
    end

    # y > 0
    function g(X)
        return [
            X[2]
        ]
    end

    i = [2]     # lower controls only x2


    # Just two variables x1 and x2
    RidgeProblem = MultiLevelProblem(2)

    RidgeProblem.addLevel!(Problem(F, G, I))
    RidgeProblem.addLevel!(Problem(f, g, i))

    println(RidgeProblem.levels())
    println(RidgeProblem.Jf(1))
    println(RidgeProblem.Jf(2))
    println(RidgeProblem[2])
    export RidgeProblem
end