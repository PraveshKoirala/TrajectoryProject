module ALEx
    include("Problems.jl")
    using .Problems

    # objective for upper
    function F(X)
        -(X[1]+X[2]+2X[3]+X[4])
    end

    # objective for lower
    function f1(X)
        -(3X[2]-X[3]-X[4]-X[1])
    end

    function f2(X)
        -(3X[3]-X[2]-X[4]-X[1])
    end

    function f3(X)
        -(3X[4]-X[2]-X[3]-X[1])
    end

    function G(X)
        []
    end

    # constraints for lower.
    function g(X)
        -[
            3X[1]+3X[2]-30,
            2X[1]+X[2]-20,
            X[3]-10,
            X[3]+X[4]-15,
            X[4]-10,
            X[1]+2X[2]+2X[3]+X[4]-40,
            -X[1], -X[2], -X[3], -X[4]
        ]
    end

    I = [1]
    i = [2,3,4]


    # Four variables
    ALProblem = MultiLevelProblem(4)

    ALProblem.addLevel!(Problem(F, G, I, 50))
    ALProblem.addLevel!(Problem((f1, f2, f3), g, i, 500))
    # add a dummy level. This won't have any effect at all and will be ignored.
    ALProblem.addLevel!(Problem(f1, g, i, 5))
    ALProblem.x_s = [0, 0, 0, 0]
    ALProblem.alpha = 1
    export ALProblem
end