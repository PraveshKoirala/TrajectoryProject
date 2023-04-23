using Ipopt
using JuMP

m = Model(Ipopt.Optimizer)

@variable(m, x)
@variable(m, y>=0)
@variable(m, λ>=0)
@constraint(m, y>=0)

w1, w2 = -5, 5.9
@objective(m, Min, (x-w1)^2 + (y-w2)^2)

δ(x) = exp(-1e12*x^2)

ϵ = 1e-6
@NLconstraint(m, -ϵ <= y-x-λ*δ(y) <= ϵ)

optimize!(m)
yx = value(x)
yy = value(y)
l = value(λ)

println(yx)
println(yy)
println(l)