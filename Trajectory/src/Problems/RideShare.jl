module RS
    include("Problems.jl")
    using .Problems
    using Plots
    using JuMP
    using Ipopt
    
    arch = "(LU)DP"
    
    # Work with Z = 3 i.e. 3 zones (Vanderbilt, Downtown, Airport)
    Z = 3
    Πᴼ = [0.2, 0.3, 0.5]                        # origin distribution
    Πᵀ = [0.05, 0.6, 0.35]                      # destination distribution
    Πᴰ = [0.1, 0.35, 0.55]                      # driver's distribution
    
    C = [ 0. 2. 8.                              # Distance matrix
          2. 0. 5.
          8. 5. 0.]
    
    OD = similar(C) .* 0
    for i in 1:Z
        except_i = setdiff(collect(1:Z), [i])
        OD[i, except_i] = Πᴼ[i] .* Πᵀ[except_i] ./ sum(Πᵀ[except_i])
    end
    OD .*= 100      # rescaled so that it works better with learning rate.
    Πᴰ .*= 10

    @assert (sum(Πᴼ) == 1 && sum(Πᵀ) == 1 && sum(Πᴰ) ≈ 10 && sum(OD) ≈ 100)

    # Dimensions
    fᵁrᵁcᵁ_dim = 1:3
    fᴸrᴸcᴸ_dim = 4:6
    Aᵁ_dim = 7:(7+Z-1)
    Aᴸ_dim = (7+Z) : (7+2Z-1)
    Pᵁ_dim = collect((7+2Z):(7+2Z+Z^2-1))
    Pᴸ_dim = collect((7+2Z+Z^2):(7+2Z+2Z^2-1))
    dim = (7+2Z+2Z^2)   # Total dimensions

    if Z != 3
        println("Don't forget to ammend the snippet below")
    end
    # Diagonal elements in demand matrices are always 0 (and unchanged)
    Pᵁ_dim_p = setdiff(Pᵁ_dim, [13, 17, 21])
    Pᴸ_dim_p = setdiff(Pᴸ_dim, [22, 26, 30])
    
    Aᴾ = [10, 10, 10]  # availability of public transit
    λ = 0.5
    g = 1               # gas cost
    fᴾ = 15             # fixed initial cost for public vehicles
    rᴾ = 0.5            # rates for public vehicles

    function F_P(X)         # Passenger's objective
        ϵ = 1e-6
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        for i in 1:Z
            Pᵁ[i, i] = 0.
            Pᴸ[i, i] = 0.
        end

        Pᴾ = reshape(OD - Pᵁ - Pᴸ, Z, Z)
        fᵁ, rᵁ, _ = X[fᵁrᵁcᵁ_dim]
        fᴸ, rᴸ, _ = X[fᴸrᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Pᵁᵢ = sum(Pᵁ, dims=2)
        Pᴸᵢ = sum(Pᴸ, dims=2)
        Pᴾᵢ = sum(Pᴾ, dims=2)
        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᵁ[i,j] * (fᵁ + rᵁ * C[i, j] + λ * Pᵁᵢ[i]/(Aᵁ[i]+ϵ)) +
                       Pᴸ[i,j] * (fᴸ + rᴸ * C[i, j] + λ * Pᴸᵢ[i]/(Aᴸ[i]+ϵ)) +
                       Pᴾ[i,j] * (fᴾ + rᴾ * C[i, j] + λ * Pᴾᵢ[i]/(Aᴾ[i]+ϵ))
            end
        end
        return acc
    end

    function G_P(X)
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴾ = reshape(OD - Pᵁ - Pᴸ, Z, Z)
        return vcat(vec(Pᵁ), vec(Pᴸ), vec(Pᴾ))
    end

    function F_D(X)         # Driver's objective
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        # Pᴾ = reshape(OD - Pᵁ - Pᴸ, Z, Z)
        _, _, cᵁ = X[fᵁrᵁcᵁ_dim]
        _, _, cᴸ = X[fᴸrᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]

        acc = 0
        for i in 1:Z
            acc += Aᵁ[i] * (cᵁ - g) * sum(Pᵁ[i, j]*C[i, j] for j in 1:Z) + 
                    Aᴸ[i] * (cᴸ - g) * sum(Pᴸ[i,j]*C[i, j] for j in 1:Z)
        end
        return -acc     # maximization
    end

    function G_D(X)
        _, _, cᵁ = X[fᵁrᵁcᵁ_dim]
        _, _, cᴸ = X[fᴸrᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        return vcat(Πᴰ - Aᵁ - Aᴸ, Aᵁ, Aᴸ, [cᵁ-g, cᴸ-g])
    end

    function F_U(X)         # Uber
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        fᵁ, rᵁ, cᵁ = X[fᵁrᵁcᵁ_dim]
        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᵁ[i,j] * (fᵁ + (rᵁ - cᵁ) * C[i, j])
            end
        end
        return -acc         # maximization
    end

    function G_U(X)
        fᵁ, rᵁ, cᵁ = X[fᵁrᵁcᵁ_dim]
        return [
            cᵁ,
            cᵁ-g
        ]
    end

    function F_L(X)         # Lyft
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        fᴸ, rᴸ, cᴸ = X[fᴸrᴸcᴸ_dim]
        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᴸ[i,j] * (fᴸ + (rᴸ - cᴸ) * C[i, j])
            end
        end
        return -acc     # maximization
    end

    function G_L(X)
        fᴸ, rᴸ, cᴸ = X[fᴸrᴸcᴸ_dim]
        return [
            cᴸ,
            cᴸ-g
        ]
    end

    function G_LU(X)
        return vcat(G_U(X), G_L(X))
    end


    function F_dummy(X)
    end
    function G_dummy(X)
    end

    # four variables
    RSProblem = MultiLevelProblem(dim)
    n = 10
    if arch == "ULDP"
        RSProblem.addLevel!(Problem(F_U, G_U, fᵁrᵁcᵁ_dim, 2n))
        RSProblem.addLevel!(Problem(F_L, G_L, fᴸrᴸcᴸ_dim, 2n))
        RSProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim), 4n))
        RSProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p), 40n))
        RSProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    elseif arch == "LUDP"
        RSProblem.addLevel!(Problem(F_L, G_L, fᴸrᴸcᴸ_dim, 2n))
        RSProblem.addLevel!(Problem(F_U, G_U, fᵁrᵁcᵁ_dim, 2n))
        RSProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim), 4n))
        RSProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p), 40n))
        RSProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    elseif arch == "(LU)DP"
        RSProblem.addLevel!(Problem((F_L, F_U), G_LU, vcat(fᴸrᴸcᴸ_dim, fᵁrᵁcᵁ_dim), 4n))
        RSProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim), 4n))
        RSProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p), 40n))
        RSProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    end

    function RSProblem.visualize(X; kwargs...)
        fᵁrᵁcᵁ = X[fᵁrᵁcᵁ_dim]
        fᴸrᴸcᴸ = X[fᴸrᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim] ./ 10
        Aᴸ = X[Aᴸ_dim] ./ 10
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z) ./ 100
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z) ./ 100
        Pᴾ = reshape(OD - Pᵁ - Pᴸ, Z, Z) ./ 100
        println("fᵁrᵁcᵁ = ", fᵁrᵁcᵁ)
        println("fᴸrᴸcᴸ = ", fᴸrᴸcᴸ)
        println("Aᵁ = ", Aᵁ)
        println("Aᴸ = ", Aᴸ)
        println("Pᵁ = ", Pᵁ)
        println("Pᴸ = ", Pᴸ)
        println("Pᴾ = ", Pᴾ)
    end

    function find_feasible_point()
        fᵁrᵁcᵁ = [5, 2, 1]
        fᴸrᴸcᴸ = [4, 1.5, 1]
        Aᵁ = Πᴰ ./ 2
        Aᴸ = Πᴰ ./ 2
        Pᵁ = OD ./ 3
        Pᴸ = OD ./ 3
        return vcat(fᵁrᵁcᵁ, fᴸrᴸcᴸ, Aᵁ, Aᴸ, vec(Pᵁ), vec(Pᴸ))
    end

    RSProblem.x_s = find_feasible_point()
    # RSProblem.x_s = [0.921729805025922, 26.205583958240894, 1.2143599111225507, 4.76329432387886, 25.011495450687793, 1.6877590743896507, 0.4008928153410466, 2.972331763578343, 0.6257519074322478, 0.034245047698801545, 0.2991469250276102, 4.799831936264837, 0.0, 0.10800499029668909, 1.2384515946854138, 0.5273070315738075, 0.0, 1.4099191544988732, 0.7537198708446567, 11.302535013863094, 0.0, 0.0, 1.0395392623228243, 2.4575180050456713, 0.07643642149219798, 0.0, 15.700322251176022, 0.1071785821558191, 0.21377907868414503, 0.0]
    RSProblem.MAX_ITER = 150
    RSProblem.alpha = 1
    export RSProblem
end