#10.1016/S0165-0114(02)00362-7
#Reported optimum f1 = 16.25 at (2.25, 0, 0, 0.25)
# After 100 iterations, obtained optimum at around -16.14585910644546 at [2.205368149901462, 0.05950914346493667, -1.0001904761908105e-8, 0.26487729336639865]

module SinhaEx1
    include("Problems.jl")
    using .Problems

    F1(X) = -(7*X[1]+3*X[2]-4*X[3]+2*X[4])      #-ve because maximization problem
    F2(X) = -(X[2]+3*X[3]+4*X[4])
    F3(X) = -(2*X[1]+X[2]+X[3]+X[4])
    G1(X) = []
    G2(X) = []
    G3(X) = -[
            X[1]+X[2]+X[3]+X[4]-5,
            X[1]+X[2]-X[3]-X[4]-2,
            -X[1]-X[2]-X[3]+1,
            -X[1]+X[2]+X[3]-1,
            X[1]-X[2]+X[3]+2*X[4]-4,
            X[1]+2*X[3]+3*X[4]-3,
            X[4]-2,
            -X[1],
            -X[2],
            -X[3],
            -X[4]
    ]
    I1 = [1, 2]      # upper controls both x1 and x2
    I2 = [3]
    I3 = [4]     # lower controls only x2


    # Four variables
    SinhaProblem = MultiLevelProblem(4)

    SinhaProblem.addLevel!(Problem(F1, G1, I1, 6))
    SinhaProblem.addLevel!(Problem(F2, G2, I2, 3))
    SinhaProblem.addLevel!(Problem(F3, G3, I3, 3))

    SinhaProblem.x_s = [0.4, 0.4, 0.4, 0.4]
    SinhaProblem.alpha = 1
    export SinhaProblem
end