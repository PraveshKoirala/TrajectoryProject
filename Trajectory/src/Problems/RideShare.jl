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

    @assert (sum(Πᴼ) == 1 && sum(Πᵀ) == 1 && sum(Πᴰ) ≈ 1 && sum(OD) ≈ 1)

    # Dimensions
    rᵁcᵁ_dim = 1:2
    rᴸcᴸ_dim = 3:4
    Aᵁ_dim = 5:(5+Z-1)
    Aᴸ_dim = (5+Z) : (5+2Z-1)
    Pᵁ_dim = collect((5+2Z):(5+2Z+Z^2-1))
    Pᴸ_dim = collect((5+2Z+Z^2):(5+2Z+2Z^2-1))
    dim = (5+2Z+2Z^2-1)   # Total dimensions

    println("Dimensions:")
    println("rᵁcᵁ_dim = ", rᵁcᵁ_dim, " (", length(rᵁcᵁ_dim), ")")
    println("rᴸcᴸ_dim = ", rᴸcᴸ_dim, " (", length(rᴸcᴸ_dim), ")")
    println("Aᵁ_dim = ", Aᵁ_dim, " (", length(Aᵁ_dim), ")")
    println("Aᴸ_dim = ", Aᴸ_dim, " (", length(Aᴸ_dim), ")")
    println("Pᵁ_dim = ", Pᵁ_dim, " (", length(Pᵁ_dim), ")")
    println("Pᴸ_dim = ", Pᴸ_dim, " (", length(Pᴸ_dim), ")")
    println()

    # Diagonal elements in demand matrices are always 0 (and unchanged)
    Pᵁ_dim_p = setdiff(Pᵁ_dim, Pᵁ_dim[1:Z+1:end])
    Pᴸ_dim_p = setdiff(Pᴸ_dim, Pᴸ_dim[1:Z+1:end])
    
    Aᴾ = [1, 1, 1]  # availability of public transit
    λ = 0.5
    g = 1               # gas cost
    fᴾ = 10             # fixed initial cost for public vehicles
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
        rᵁ, _ = X[rᵁcᵁ_dim]
        rᴸ, _ = X[rᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Pᵁᵢ = sum(Pᵁ, dims=2)
        Pᴸᵢ = sum(Pᴸ, dims=2)
        Pᴾᵢ = sum(Pᴾ, dims=2)
        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᵁ[i,j] * (rᵁ * C[i, j] + λ * Pᵁᵢ[i]/(Aᵁ[i]+ϵ)) +
                       Pᴸ[i,j] * (rᴸ * C[i, j] + λ * Pᴸᵢ[i]/(Aᴸ[i]+ϵ)) +
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
        _, cᵁ = X[rᵁcᵁ_dim]
        _, cᴸ = X[rᴸcᴸ_dim]
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
        _, cᵁ = X[rᵁcᵁ_dim]
        _, cᴸ = X[rᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        return vcat(Πᴰ - Aᵁ - Aᴸ, Aᵁ, Aᴸ, [cᵁ-g, cᴸ-g])
    end

    function F_U(X)         # Uber
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        rᵁ, cᵁ = X[rᵁcᵁ_dim]
        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᵁ[i,j] * ((rᵁ - cᵁ) * C[i, j])
            end
        end
        return -acc         # maximization
    end

    function G_U(X)
        rᵁ, cᵁ = X[rᵁcᵁ_dim]
        return [
            cᵁ,
            cᵁ-g
        ]
    end

    function F_L(X)         # Lyft
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        rᴸ, cᴸ = X[rᴸcᴸ_dim]
        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᴸ[i,j] * ((rᴸ - cᴸ) * C[i, j])
            end
        end
        return -acc     # maximization
    end

    function G_L(X)
        rᴸ, cᴸ = X[rᴸcᴸ_dim]
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
        RSProblem.addLevel!(Problem(F_U, G_U, rᵁcᵁ_dim, 2n))
        RSProblem.addLevel!(Problem(F_L, G_L, rᴸcᴸ_dim, 2n))
        RSProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim), 4n))
        RSProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p), 40n))
        RSProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    elseif arch == "LUDP"
        RSProblem.addLevel!(Problem(F_L, G_L, rᴸcᴸ_dim, 2n))
        RSProblem.addLevel!(Problem(F_U, G_U, rᵁcᵁ_dim, 2n))
        RSProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim), 4n))
        RSProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p), 40n))
        RSProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    elseif arch == "(LU)DP"
        RSProblem.addLevel!(Problem((F_L, F_U), G_LU, vcat(rᴸcᴸ_dim, rᵁcᵁ_dim), 4n))
        RSProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim), 4n))
        RSProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p), 40n))
        RSProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    end

    function RSProblem.visualize(X; kwargs...)
        rᵁcᵁ = X[rᵁcᵁ_dim]
        rᴸcᴸ = X[rᴸcᴸ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴾ = reshape(OD - Pᵁ - Pᴸ, Z, Z)
        println("rᵁcᵁ = ", rᵁcᵁ)
        println("rᴸcᴸ = ", rᴸcᴸ)
        println("Aᵁ = ", Aᵁ)
        println("Aᴸ = ", Aᴸ)
        println("Pᵁ = ", Pᵁ)
        println("Pᴸ = ", Pᴸ)
        println("Pᴾ = ", Pᴾ)
    end

    function find_feasible_point()
        rᵁcᵁ = [2, 1]
        rᴸcᴸ = [1.5, 1]
        Aᵁ = Πᴰ ./ 2
        Aᴸ = Πᴰ ./ 2
        Pᵁ = OD ./ 3
        Pᴸ = OD ./ 3
        return vcat(rᵁcᵁ, rᴸcᴸ, Aᵁ, Aᴸ, vec(Pᵁ), vec(Pᴸ))
    end

    RSProblem.x_s = find_feasible_point()
    # RSProblem.x_s = [0.921729805025922, 26.205583958240894, 1.2143599111225507, 4.76329432387886, 25.011495450687793, 1.6877590743896507, 0.4008928153410466, 2.972331763578343, 0.6257519074322478, 0.034245047698801545, 0.2991469250276102, 4.799831936264837, 0.0, 0.10800499029668909, 1.2384515946854138, 0.5273070315738075, 0.0, 1.4099191544988732, 0.7537198708446567, 11.302535013863094, 0.0, 0.0, 1.0395392623228243, 2.4575180050456713, 0.07643642149219798, 0.0, 15.700322251176022, 0.1071785821558191, 0.21377907868414503, 0.0]
    RSProblem.MAX_ITER = 150
    RSProblem.alpha = fill(1., dim)
    RSProblem.alpha[Aᵁ_dim] .= 0.1
    RSProblem.alpha[Aᴸ_dim] .= 0.1
    RSProblem.alpha[Pᵁ_dim] .= 0.01
    RSProblem.alpha[Pᴸ_dim] .= 0.01
    export RSProblem
end