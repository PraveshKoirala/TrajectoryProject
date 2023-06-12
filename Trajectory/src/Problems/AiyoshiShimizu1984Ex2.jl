module AS1984
    include("Problems.jl")
    using .Problems

    function F(x)
        2x[1] + 2x[2] - 3x[3] - 3x[4] -60
    end

    # Unconstrained
    function G(x)
        -[
            x[1]+x[2]+x[3]-2x[4]-40, 
            x[1]-50,
            x[2]-50,
            -x[1],
            -x[2]
        ]
    end
    
    I = [1, 2]      # upper controls both x1 and x2

    function f(x)
        (x[3]-x[1]+20)^2 + (x[4]-x[2]+20)^2
    end

    # y > 0
    function g(x)
        -[
            2x[3]-x[1]+10, 
            2x[4]-x[2]+10, 
            -x[3]-10,-x[4]-10, 
            x[3]-20, 
            x[4]-20
        ]
    end

    i = [3, 4]     # lower controls only x2


    AS1984Problem = MultiLevelProblem(4)

    AS1984Problem.addLevel!(Problem(F, G, I))
    AS1984Problem.addLevel!(Problem(f, g, i))

    AS1984Problem.x_s = [10, 10, 20, 20]
    export AS1984Problem
end