using Ipopt
# include("Problems/Sinha.jl")
include("Problems/Tilahun.jl")
# using .SinhaEx1: SinhaProblem as p
using .Tilahun: TilahunProblem as p

# include("Problems/NestedToll.jl")

using JuMP
import Symbolics
using Ipopt

function solvefull(p, x_s, η=1)
    """
    p is the multilevel problem and x_s is the tentative solution. We now iteratively construct
    """
    ϵ=1e-6
    dim = length(p.X)
    
    A = []
    In = []
    I = vcat([], p[p.levels()].I)
    duals = [1]
    for j = p.levels()-1 : -1 : 1
        # x is a symbolic variable
        @Symbolics.variables x[1:dim]
        O = p[j+1].f(x)[end]            # players objectives and constraints
        C = p[j+1].g(x)
        values = p[j+1].g(x_s)
        if any(p[j+1].g(x_s) .< -ϵ) throw(DomainError("none of the constraints can be negative")) end
        A = vcat(A, C[-ϵ .<= values .<= ϵ])                     # active constraints (close to 0)
        In= vcat(In, C[values .>= ϵ] .- ϵ)                             # inactive constraints (enforce strict > by adding ϵ)
        # Implicitly embed bottom constraints
        if length(A) == 0
            # so there are constraints but none of them are active at this point. Which means x_s should be at the interior
            # that means that objective must be optimal and the optimality must be maintained.
            In = push!(In, p[j+1].f(x_s)[end]-O)       # this must be == 0
        end

        # Add all top constraints
        O = p[j].f(x)[end]            # players objectives and constraints
        C = p[j].g(x)
        In = vcat(In, p[j].g(x))

        # x stops being a symbolic variable and becomes a JuMP variable. Super confusing.
        global m = Nothing
        global x = Nothing
        m = Model(Ipopt.Optimizer)
        # set_silent(m)
        @variable(m, x[i=1:dim], start=x_s[i])
        I = vcat(I, p[j].I)                                # variables that this player control
        Ī = setdiff(1:dim, I)                   # player controls their and all lower variables
        @constraint(m, x[Ī] .== x_s[Ī])
        @constraint(m, x_s.-η .<= x .<= x_s.+η)    # specify bounds
        
        # Add constraints
        for ci in eachindex(A)
            eval(Meta.parse("@NLconstraint(m, $(repr(A[ci])) == 0)"))  
        end

        for ci in eachindex(In)
            eval(Meta.parse("@NLconstraint(m, $(repr(In[ci])) >= 0)"))
        end
        eval(Meta.parse("@NLobjective(m, Min, $(repr(O)))" ))    # Solve

        optimize!(m)
        x_s = value.(x)
        # for con in all_constraints(m; include_variable_in_set_constraints = true)
        #     duals dual(con))
        # end
        println(x_s)
    end
    return x_s
end

# solvefull(p, [2.18, 0.09, -1.000e-8, 0.27])
solvefull(p, [0.44, 0, 0.44])       # for Tilahun