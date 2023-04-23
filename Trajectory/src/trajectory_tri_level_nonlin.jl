using Symbolics

# Setup the environment
r = 5        #radius of the feasibility region in positive quadrant
R = [r, r]  # center of the feasibility circle

O = [12, 9]  # obstacle center
or = 2      # radius of the obstacle

N_τ = 20     #number of trajectory points is n+1
D = 20      #destination plane, any value larger than the obstacle should be fine

function NLinearPI(xnm1)
    δ = 1
    # Linear policy, just move δ step forward in x direction
    return xnm1[1] + δ, xnm1[2] + 0.5*sin(xnm1[1]+δ)
end

PI = NLinearPI

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

F = [F_tau, F_t, F_T]
G = vcat(G_t(x_), G_T(x_))

# Level to decision variable mapping
decision_variables = Dict(1 => xdim, 2 => vcat(xdim, tdim), 3 => Tdim)

# Level to constraint mapping
G_t_length = length(G_t(x_))
G_T_length = length(G_T(x_))
G_length = G_t_length + G_T_length
constraints_map = Dict(1=>[], 2=>1:G_t_length, 3=>(G_t_length+1):(G_t_length+G_T_length))

# Manual initial feasible point.
x_s = [5., 0., 100., 100.];