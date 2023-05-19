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

    function f(X)
        return 0.5 * (X[1]-X[2])^2
    end

    # y > 0
    function g(X)
        return [
            X[2]
        ]
    end


    # Just two variables x1 and x2
    RidgeProblem = MultiLevelProblem(2)

    RidgeProblem.addLevel!(Problem(F, G))
    RidgeProblem.addLevel!(Problem(f, g))

    println(RidgeProblem.levels())
    println(RidgeProblem.Jf(1))
    println(RidgeProblem.Jf(2))
    println(RidgeProblem[2])
    export RidgeProblem
end