using Symbolics

# Setup the environment
r = 5        #radius of the feasibility region in positive quadrant
R = [r, r]  # center of the feasibility circle

O = [12, 9]  # obstacle center
or = 2      # radius of the obstacle

N_τ = 20     #number of trajectory points is n+1
D = 20      #destination plane, any value larger than the obstacle should be fine

function LinearPI(xnm1)
    δ = 1
    # Linear policy, just move δ step forward in x direction
    return xnm1[1] + δ, xnm1[2]
end

PI = LinearPI

# Tri-level problem
L = 3

dim = 4;
xdim = 1:2;
tdim = 3:3;
Tdim = 4:4;

# x is upper level, y is lower level. _ means that they are symbolic variables
@Symbolics.variables x_[1:dim]
x_ = Symbolics.scalarize(x_)

# Slack from the obstacle. >=0 if trajectory is not touching the obstacle.
function g_tau(tau)
    return (tau[1]-O[1])^2 + (tau[2]-O[2])^2 - or^2
end

function chi(x0)
    # check if the initial point in feasible region
    return r^2-(x0[1]-R[1])^2-(x0[2]-R[2])^2
end

function F_T(x)
    return -x[Tdim][1]
end

function G_T(x)
    T = x[Tdim][1]
    constraints = [g_tau(x) - T];
    x_next = x;
    # Generate next trajectory points using the policy
    for _ = 2:N_τ
        x_next = PI(x_next)
        push!(constraints, g_tau(x_next) - T)
    end
    return constraints
end

# objective for upper
function F_t(x)
    t = x[tdim][1]
    return 1/2 * t^2
end

function G_t(x)
    T = x[Tdim][1];
    t = x[tdim][1]
    return [
        t-T,
        chi(x),
        T
    ]
end


function F_tau(x)
    return - (D - x[1])
end


using JuMP
using Ipopt

function find_feasible_point(x_s)
    # Use Ipopt to find an initial feasible point.
    # The point will be inside the constraint boundary but not exactly feasible in Multilevel optimization sense
    # After a point is obtained, we must let the lower level players optimize to truly obtain a feasible starting point.
    
    # To find a starting point, we linearly weigh the function objectives with the last level getting the most
    # weight. This is in hopes that we'll achieve something sensible or near the local optimum.
    r1 = 1;
    r2 = 10;
    r3 = 100;

    m = Model(Ipopt.Optimizer)
    @variable(m, x[i=1:dim], start=x_s[i])

    C = vcat(G_T(x), G_t(x))
    lc = length(C)

    @variable(m, s[1:lc] .>= 0)
    @NLconstraint(m, [i=1:lc], C[i]>=s[i])

    obj = @expression(m, r1 * F_tau(x) + r2 * F_t(x) + r3 * F_T(x));
    @NLobjective(m, Min, obj)
    optimize!(m)
    return (termination_status(m), value.(x))
end

# Find a feasible point. Start with some random values
x_s = [0, 0, 1, 2];

status, x_s = find_feasible_point(x_s)

x_s = [0.015401027606127604, 5.14272294205888, 11, 10.878760000000172];

println(status)
println(x_s)