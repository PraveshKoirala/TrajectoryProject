module CalamaiVicente1994b
    include("Problems.jl")
    using .Problems

    # objective for upper
    function F(x)
        y = x[5:end]
        ay = y[1]+y[3]  # a.y
        by = y[2]+y[4]  #
        -(200 .- ay).*(ay)-(160 .- by).*(by)
    end

    # objective for lower
    function f(x)
        y = x[5:end]
        return (y[1]-4)^2 + (y[2]-13)^2 + (y[3]-35)^2 + (y[4]-2)^2
    end

    # constraints for upper BOLIB has the form g(x,y)<=0 which is annoying, so we'll return -G to make it >= 0
    function G(x)
        -vec([sum(x[1:4])-40; x[1:4]-[10 5 15 20]';-x[1:4] ])
    end

    # constraints for lower.
    function g(x)
        y = x[5:end]
        A = [0.4 0.7 0 0; 0.6 0.3 0 0;0 0 0.4 0.7; 0 0 0.6 0.3];
        -vec([A*y-x[1:4]; y-[20 20 40 40]'; -y])
    end

    I = [1, 2, 3, 4]
    i = [5, 6, 7, 8]


    Bard1988Ex2Problem = MultiLevelProblem(8)

    Bard1988Ex2Problem.addLevel!(Problem(F, G, I))
    Bard1988Ex2Problem.addLevel!(Problem(f, g, i))

    Bard1988Ex2Problem.x_s = [0, 0, 0, 0, 0, 0, 0, 0];
    export Bard1988Ex2Problem
end