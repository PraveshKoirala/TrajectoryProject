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
    end

    (F::FunctionTuple)(x) = Tuple(f(x) for f in F)
    Problem(f::Function, g, I) = Problem((f,), g, I)

    @kwdef mutable struct MultiLevelProblem
        N::Int  # Number of variables
        X::Vector{Num}
        P::Vector{Problem} = Vector{Problem}()
        levels::Function = @memoize L() -> length(P)

        # Jacobian of the objective function for i-th level
        Jf = @memoize J1(i) -> Symbolics.gradient(P[i].f(X), X)
        # Jacobian of the constraint functions for the i-th level
        Jg = @memoize J2(i) -> Symbolics.jacobian(P[i].g(X), X)
        
        addLevel! = (p::Problem) -> push!(P, p);

        # x_s is initial feasible point. Default is all zeros
        x_s :: Vector{Float64} = zeros(N)
        visualize:: Function = (x; kwargs...) -> ()   # the visualization function

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