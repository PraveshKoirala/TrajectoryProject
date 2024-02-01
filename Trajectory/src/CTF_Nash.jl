module CTF
    include("constants.jl")
    include("Problems.jl")
    include("spline_methods.jl")
    using .Problems

    xmin = -5
    xmax = 5
    ymin = -5
    ymax = 5

    
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
                                    longlength, true)
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
                                    longlength, true)
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
    CTFProblem = MultiLevelProblem(2)

    CTFProblem.addLevel!(Problem((F1, F2), G1, vcat(I1, I2), 500, 4))
    CTFProblem.addLevel!(Problem((F3,), G3, I3, 0, 0))    # Dummy level, does not do anything.

    CTFProblem.x_s = rand(4*K) .* (xmax-xmin) .+ xmin

    
    # CTFProblem.x_s[1] = start1[1]
    # CTFProblem.x_s[2] = start1[2]

    # CTFProblem.x_s[2*K+1] = start2[1]
    # CTFProblem.x_s[2*K+2] = start2[2]

    CTFProblem.x_s = vcat(repeat(start1, K), repeat(start2, K))
    CTFProblem.MAX_ITER = 100
    export CTFProblem
end