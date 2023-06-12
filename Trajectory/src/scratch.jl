using JuMP, Ipopt, HiGHS


# todo: use the start to start closer
w = [-8, 1]
W = [1e6, 1e-6]
M = 120
ϕ = 10

x_s = [0, 0]
for i = 1:12
    global x_s
    # Solve full for first and last tries
    if i == 1 || i == 12
        MAX_ITER = 20
    else
        MAX_ITER = 2
    end

    m = Model(HiGHS.Optimizer)
    set_optimizer_attribute(m, "max_iter", MAX_ITER)
    set_silent(m)
    @variable(m, x[i=1:2], start=x_s[i])
    f1 = @expression(m, (x[1]-w[1])^2 + (x[2]-w[2])^2)
    f2 = @expression(m, (x[1]-x[2])^2)
    @constraint(m, x[2]>=0)
    @objective(m, Min, W[1]*f1 + W[2]*f2)
    W[1] /= ϕ
    W[2] *= ϕ
    optimize!(m)
    x_s = value.(x)
    println("Obtained value", x_s)
end
