module NToll
    include("Problems.jl")
    using .Problems
    # X = [t1, t2, p1, p2, p3]
    # D = 6       # obtain [2.244008132986424, 4.0773767557423515, 0.03939813363997559, 0.960595919057888, 5.9377111648188955e-6]
                # True answer [>=2, 4, 0, 1, 0]
    D = -1.5  # obtains [1.2898016970530486, 0.6050192328926055, 0.21653359688861815, -5.294905701298197e-9, 0.7834664070865699]
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

    NTollProblem.addLevel!(Problem(F1, G1, I1, 5))
    NTollProblem.addLevel!(Problem(F2, G2, I2, 5))
    NTollProblem.addLevel!(Problem(F3, G3, I3, 5))

    NTollProblem.x_s = [0, 0, 1, 0, 0]
    NTollProblem.alpha = 0.2
    export NTollProblem
end