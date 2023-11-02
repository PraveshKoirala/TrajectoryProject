module NRidge
    include("Problems.jl")
    using .Problems
    using Random
    Random.seed!(20)
    # Number of levels
    N = 5
    
    NS = 6

    w = rand(0:10, N)
    println(w)
    # Just two variables x1 and x2
    NRidgeProblem = MultiLevelProblem(N)
    NRidgeProblem.MAX_ITER = 25
    NRidgeProblem.x_s = zeros(N)

    norm(v) = sum(v[i]^2 for i in eachindex(v))
    
    function get_G(ix)
        G = (X)->[]
        if ix>1
            G = (X)->[X[ix-1]-X[ix]]
        end
        return G
    end

    function get_F(n)
        F(X) = norm(X[n:end] - w[n:end])
    end
    for n = 1:N
        NRidgeProblem.addLevel!(Problem(get_F(n), get_G(n), [n], NS))
    end
    
    NRidgeProblem.alpha = 5
    function NRidgeProblem.visualize(x; kwargs...)
        println(x)
        return;
        println("Initialized with value of w = ", w)
        # ideal answer
        mean(X) = sum(X)/length(X)
        function ideal(w, i, e0)
            if i == N return min(w[i], e0) end
            ans = min(e0, mean(w[i:end]))
            return vcat([ans], ideal(w, i+1, ans))
        end
        function error(X)
            return norm(w-X)
        end
        xstar = ideal(w, 1, 9999)
        println("Predicted Ideal X* = ", xstar)
        println("Ideal Error = ", error(xstar))
        println("Obtained Error = ", error(x))
    end
    export NRidgeProblem
end