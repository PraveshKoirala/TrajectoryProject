"""
A module that provides an interface for a multi-level problem.
"""
module Problems
    using Symbolics
    using Memoize
    import Base

    const FunctionTuple{N} = NTuple{N, Function} where N
    struct Problem
        f::FunctionTuple    # objective function
        g::Function         # the Constraint function
        I::Vector{UInt}     # Index set of variables
        n::Int              # Number of samples to generate
    end

    (F::FunctionTuple)(x) = Tuple(f(x) for f in F)
    Problem(f::Function, g, I, n) = Problem((f,), g, I, n)

    @kwdef mutable struct MultiLevelProblem
        N::Int  # Number of variables
        X::Vector{Num}
        P::Vector{Problem} = Vector{Problem}()
        levels::Function = () -> length(P)

        # Jacobian of the objective function for i-th level
        Jf = (i) -> Symbolics.gradient(P[i].f(X), X)
        # Jacobian of the constraint functions for the i-th level
        Jg = (i) -> Symbolics.jacobian(P[i].g(X), X)
        
        addLevel! = (p::Problem) -> push!(P, p);

        # x_s is initial feasible point. Default is all zeros
        x_s :: Vector{Float64} = zeros(N)
        visualize:: Function = (x; kwargs...) -> ()   # the visualization function
        alpha = 2
        cooldown = 1        # means that the alpha will remain constant throughout
        MAX_ITER = 100      # default iteration of 100
    end

    Base.getindex(MLP::MultiLevelProblem, i) = MLP.P[i]
    

    function MultiLevelProblem(N::Int)
        @variables X[1:N]
        P = MultiLevelProblem(
            N = N, 
            X = collect(X)
        )
        return P
    end
    
    export Problem, MultiLevelProblem
end