"""
A module that provides an interface for a multi-level problem.
"""
module Problems
    using Symbolics
    using Memoize
    import Base

    struct Problem
        f::Function         # the objective function
        g::Function         # the Constraint function
    end

    function evaluate(expression, x)
        vars_mapping = Dict(zip(x_, x))
        return substitute(expression, vars_mapping)
    end

    @kwdef struct MultiLevelProblem
        N::Int  # Number of levels
        X::Vector{Num}
        P::Vector{Problem} = Vector{Problem}()
        levels::Function = () -> length(P)

        # Jacobian of the objective function for i-th level
        @memoize Jf = (i) -> Symbolics.gradient(P[i].f(X), X)
        # Jacobian of the constraint functions for the i-th level
        @memoize Jg = (i) -> Symbolics.jacobian(P[i].g(X), X)

        addLevel! = (p::Problem) -> push!(P, p);
    end

    Base.getindex(MLP::MultiLevelProblem, i) = MLP.P[i]
    
    function MultiLevelProblem(N::Int)
        X = Symbolics.variables(:x, 1:N)
        P = MultiLevelProblem(
            N=N, 
            X=X
        )
        return P
    end
    
    export Problem, MultiLevelProblem
end