module CTF
    include("Problems.jl")
    include("spline_methods.jl")
    using .Problems

    xmin = -5
    xmax = 5
    ymin = -5
    ymax = 5

    flags = [
             0 2;
             2 2;]
    K = 4   # number of knot points x0, x1, x2, x3 back to x0 (not counted)

    function F1(X)
        knots1 = reshape(X[1:(2*K)], (2, K))'
        knots2 = reshape(X[2*K+1:4*K], (2, K))'
        grid1 = get_grid(knots1)
        grid2 = get_grid(knots2)
        t1 = get_coords(knots1, grid1)
        t2 = get_coords(knots2, grid2)
        l1 = get_length(t1)
        l2 = get_length(t2)
        id = l1 < l2 ? "S" : "L"
        short, shortlength, long, longlength = l1 < l2 ? (t1, l1, t2, l2) : (t2, l2, t1, l1)
        short_score, long_score = score(flags, short, shortlength, long,
                                    longlength, false, 1)
        if id == "S" return -short_score else return -long_score end
    end
    
    function F2(X)
        knots1 = reshape(X[2*K+1:4*K], (2, K))'
        knots2 = reshape(X[1:(2*K)], (2, K))'
        grid1 = get_grid(knots1)
        grid2 = get_grid(knots2)
        t1 = get_coords(knots1, grid1)
        t2 = get_coords(knots2, grid2)
        l1 = get_length(t1)
        l2 = get_length(t2)
        id = l1 < l2 ? "S" : "L"
        short, shortlength, long, longlength = l1 < l2 ? (t1, l1, t2, l2) : (t2, l2, t1, l1)
        short_score, long_score = score(flags, short, shortlength, long,
                                    longlength, true, 1)
        if id == "S" return -short_score else return -long_score end
    end

    function F3(X)
    end

    # Constraints are implicitly defined by objective function
    G1(X) = []
    G2(X) = [] 
    G3(X) = []
    
    I1 = collect(3:(2*K))      # x,y for all knots
    I2 = collect((2*K+3):4*K)
    I3 = []     

    # Three variables
    CTFProblem = MultiLevelProblem(3)

    CTFProblem.addLevel!(Problem((F1,), G1, I1, 10, 1))
    CTFProblem.addLevel!(Problem((F2,), G2, I2, 20, 2))
    CTFProblem.addLevel!(Problem((F3,), G3, I3, 0, 0))    # Dummy level, does not do anything.

    CTFProblem.x_s = rand(4*K) .* (xmax-xmin) .+ xmin

    start1 = [1, -2]
    start2 = [-1, -2]
    # CTFProblem.x_s[1] = start1[1]
    # CTFProblem.x_s[2] = start1[2]

    # CTFProblem.x_s[2*K+1] = start2[1]
    # CTFProblem.x_s[2*K+2] = start2[2]

    CTFProblem.x_s = vcat(repeat(start1, 4), repeat(start2, 4))
    CTFProblem.MAX_ITER = 300
    export CTFProblem
end