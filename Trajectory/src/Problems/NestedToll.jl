module NToll
    include("Problems.jl")
    using .Problems
    # X = [t1, t2, p1, p2, p3]
    D = 6       # obtain -4.001789931486663 at [2.387336379978474, 4.032828339454866, 0.005553292343221871, 0.9890161491937394, 0.005430548871248347]
                # True answer [>=2, 4, 0, 1, 0]
    # D = -1.5  # obtains -0.25447935481193085 at [1.2807633820062718, 0.14966693277741025, 0.19869349125269498, 4.658674784165937e-8, 0.8013064601924402]
    #           # True answer [1, >= 0, 0.25, 0, 0.75]
    F1(X) = -(X[3] * X[1] + X[4] * X[2])
    F2(X) = X[3] * (X[1]+X[3]) + (X[4]+X[5])*(X[4]+X[5])
    F3(X) = X[4] * (X[2]+X[4]) + X[5]*(X[5]+D)
    G1(X) = []
    G2(X) = []
    G3(X) = [
           X[3]+X[4]+X[5]-1,
           1-X[3]-X[4]-X[5],
           X[1], X[2], X[3], X[4], X[5]
    ]
    I1 = [1, 2] 
    I2 = [3]
    I3 = [4, 5]

    # Three variables
    NTollProblem = MultiLevelProblem(5)

    NTollProblem.addLevel!(Problem(F1, G1, I1, 7))
    NTollProblem.addLevel!(Problem(F2, G2, I2, 7))
    NTollProblem.addLevel!(Problem(F3, G3, I3, 5))

    NTollProblem.x_s = [0, 0, 1, 0, 0]
    NTollProblem.alpha = 0.15
    export NTollProblem
end