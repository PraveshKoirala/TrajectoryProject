#10.1016/S0165-0114(02)00362-7
#Reported optimum f1 = 16.25 at (2.25, 0, 0, 0.25)
# After 100 iterations, obtained optimum at around -16.14585910644546 at [2.205368149901462, 0.05950914346493667, -1.0001904761908105e-8, 0.26487729336639865]

module COExample
    include("Problems.jl")
    using .Problems

    F1(X) = (X[1]-X[2])*(X[1]-X[2])
    F2(X) = -(X[1]-X[2])*(X[1]-X[2])
    G1(X) = []
    G2(X) = [X[2], 10-X[2]]
    
    I1 = [1]      # upper controls both x1 and x2
    I2 = [2]


    # Four variables
    COProblem = MultiLevelProblem(2)

    COProblem.addLevel!(Problem(F1, G1, I1, 5))
    COProblem.addLevel!(Problem(F2, G2, I2, 5))

    COProblem.x_s = [10, 10]
    COProblem.alpha = 0.5
    export COProblem
end