# A New Algorithm for Multilevel Optimization Problems 
# Using Evolutionary Strategy, Inspired 
# by Natural Adaptation
# Reported approximate optimum f1 = -0.5 at (0.5, 0, 0.0095) (Completely wrong)
# After 100 iterations, obtained optimum at around -0.49798952448804745 at [0.49798952448804745, 1.083378936092657e-18, 0.4979895319788691]

# O (does not, wrong) -> P1 (matches my answer)
module Tilahun
    include("Problems.jl")
    using .Problems

    F1(X) = -X[1]+4*X[2]
    F2(X) = 2*X[2]+X[3]
    F3(X) = -X[3]^2+X[2]
    G1(X) = [1-X[1]-X[2]]
    G2(X) = [
        -X[3]+2X[1]-X[2]]
    G3(X) = [
            X[1]-X[3],
            X[3], 1-X[3],
            X[1], 0.5-X[1],
            X[2], 1-X[2]
    ]
    I1 = [1]      # upper controls both x1 and x2
    I2 = [2]
    I3 = [3]     # lower controls only x2

    # Three variables
    TilahunProblem = MultiLevelProblem(3)

    TilahunProblem.addLevel!(Problem(F1, G1, I1, 5))
    TilahunProblem.addLevel!(Problem(F2, G2, I2, 5))
    TilahunProblem.addLevel!(Problem(F3, G3, I3, 5))

    TilahunProblem.x_s = [0, 0, 0]
    TilahunProblem.alpha = 0.2
    export TilahunProblem
end