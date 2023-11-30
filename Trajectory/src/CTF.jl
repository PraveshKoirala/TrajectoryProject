module CTF
    include("Problems.jl")
    using .Problems

    xmin = -5
    xmax = 5
    ymin = -5
    ymax = 5

    flags = [1 1;
             2 1;
             3 2;
             3 3]
    K = length(flags)

    nflags = size(flags)[1]
    points = 100
    grid = hcat([(1 .- t).^(K-n) .* t.^(n) .* bins[n+1] for n in 0:K]...)

    function score(T1, T2)
        # calculate score for T1 w.r.t T2
        
    end
    function f1(X)
        x = X[1:knots*2]
        y = X[2*knots+1:end]

    end
    
    I1 = [1]      # upper controls both x1 and x2
    I2 = [2]
    I3 = [3]     # lower controls only x2

    # Three variables
    CTFProblem = MultiLevelProblem(3)

    CTFProblem.addLevel!(Problem(F1, G1, I1, 5))
    CTFProblem.addLevel!(Problem(F2, G2, I2, 5))
    CTFProblem.addLevel!(Problem(F3, G3, I3, 5))

    CTFProblem.x_s = [0, 0, 0]
    CTFProblem.alpha = 0.2
    export CTFProblem
end