# Use for debugging in VS Code
using Random
using Trajectory
Random.seed!(27)

t = @elapsed run_opt()
println("Total Time taken: ", t)