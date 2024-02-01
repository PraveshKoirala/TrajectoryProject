module RSI
    include("Problems.jl")
    include("RideShare.jl")
    using .Problems
    using Plots
    using JuMP
    using Ipopt
    using .RS: F_P as F_P0, F_D as F_D0
    
    arch = "(LUI)DP"
    
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
    rᴵ_dim = 5:5        # only rate is defined for indrive as commission would just be (1-γ)r
    Aᵁ_dim = 6:(6+Z-1)
    Aᴸ_dim = Aᵁ_dim .+ Z
    Aᴵ_dim = Aᴸ_dim .+ Z
    Pᵁ_dim = collect(Aᴵ_dim[end]+1 : Aᴵ_dim[end]+ Z^2)
    Pᴸ_dim = collect(Pᵁ_dim .+ Z^2)
    Pᴵ_dim = collect(Pᴸ_dim .+ Z^2)
    dim = Pᴵ_dim[end]   # Total dimensions

    println("Dimensions:")
    println("rᵁcᵁ_dim = ", rᵁcᵁ_dim, " (", length(rᵁcᵁ_dim), ")")
    println("rᴸcᴸ_dim = ", rᴸcᴸ_dim, " (", length(rᴸcᴸ_dim), ")")
    println("rᴵ_dim = ", rᴵ_dim, " (", length(rᴵ_dim), ")")
    println("Aᵁ_dim = ", Aᵁ_dim, " (", length(Aᵁ_dim), ")")
    println("Aᴸ_dim = ", Aᴸ_dim, " (", length(Aᴸ_dim), ")")
    println("Aᴵ_dim = ", Aᴵ_dim, " (", length(Aᴵ_dim), ")")
    println("Pᵁ_dim = ", Pᵁ_dim, " (", length(Pᵁ_dim), ")")
    println("Pᴸ_dim = ", Pᴸ_dim, " (", length(Pᴸ_dim), ")")
    println("Pᴵ_dim = ", Pᴵ_dim, " (", length(Pᴵ_dim), ")")
    println()

    # Diagonal elements in demand matrices are always 0 (and unchanged)
    Pᵁ_dim_p = setdiff(Pᵁ_dim, Pᵁ_dim[1:Z+1:end])
    Pᴸ_dim_p = setdiff(Pᴸ_dim, Pᴸ_dim[1:Z+1:end])
    Pᴵ_dim_p = setdiff(Pᴵ_dim, Pᴵ_dim[1:Z+1:end])

    Aᴾ = [1, 1, 1]  # availability of public transit
    λ = 0.5
    γ = 0.1             # in-drive's commission
    g = 1               # gas cost
    fᴾ = 10             # fixed initial cost for public vehicles
    rᴾ = 0.5            # rates for public vehicles

    function F_P(X)         # Passenger's objective
        ϵ = 1e-6
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z)
        Pᴾ = reshape(OD - Pᵁ - Pᴸ - Pᴵ, Z, Z)

        rᵁ, _ = X[rᵁcᵁ_dim]
        rᴸ, _ = X[rᴸcᴸ_dim]
        rᴵ = X[rᴵ_dim][1]

        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]

        Pᵁᵢ = sum(Pᵁ, dims=2)
        Pᴸᵢ = sum(Pᴸ, dims=2)
        Pᴾᵢ = sum(Pᴾ, dims=2)
        Pᴵᵢ = sum(Pᴵ, dims=2)

        acc = 0
        for i in 1:Z
            for j in 1:Z
                acc += Pᵁ[i,j] * (rᵁ * C[i, j] + λ * Pᵁᵢ[i]/(Aᵁ[i]+ϵ)) +
                       Pᴸ[i,j] * (rᴸ * C[i, j] + λ * Pᴸᵢ[i]/(Aᴸ[i]+ϵ)) +
                       Pᴵ[i,j] * (rᴵ * C[i, j] + λ * Pᴵᵢ[i]/(Aᴵ[i]+ϵ)) +
                       Pᴾ[i,j] * (fᴾ + rᴾ * C[i, j] + λ * Pᴾᵢ[i]/(Aᴾ[i]+ϵ))
            end
        end
        return acc
    end

    function G_P(X)
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z)
        Pᴾ = reshape(OD - Pᵁ - Pᴸ - Pᴵ, Z, Z)
        return vcat(vec(Pᵁ), vec(Pᴸ), vec(Pᴾ), vec(Pᴵ))
    end

    function F_D(X)         # Driver's objective
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z)
        # Pᴾ = reshape(OD - Pᵁ - Pᴸ, Z, Z)
        _, cᵁ = X[rᵁcᵁ_dim]
        _, cᴸ = X[rᴸcᴸ_dim]
        rᴵ = X[rᴵ_dim][1]
        cᴵ = (1-γ)rᴵ

        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]

        acc = 0
        for i in 1:Z
            acc += Aᵁ[i] * (cᵁ - g) * sum(Pᵁ[i, j]*C[i, j] for j in 1:Z) + 
                    Aᴸ[i] * (cᴸ - g) * sum(Pᴸ[i,j]*C[i, j] for j in 1:Z) + 
                    Aᴵ[i] * (cᴵ - g) * sum(Pᴵ[i, j]*C[i, j] for j in 1:Z)
        end
        return -acc     # maximization
    end

    function G_D(X)
        _, cᵁ = X[rᵁcᵁ_dim]
        _, cᴸ = X[rᴸcᴸ_dim]
        rᴵ = X[rᴵ_dim][1]
        cᴵ = (1-γ) * rᴵ
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]

        return vcat(Πᴰ - Aᵁ - Aᴸ - Aᴵ, Aᵁ, Aᴸ, Aᴵ, [cᵁ-g, cᴸ-g, cᴵ-g])
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

    # Indrive
    function F_I(X)
        util_D = -(1-γ) * F_D(X) # utility for driver 
        util_P = -F_P(X) # utility for passenger
        return -(util_D - base_D)*(util_P - base_P)  # max
     end
 
     function G_I(X)
         rᴵ = X[rᴵ_dim][1]
         cᴵ = (1-γ) * rᴵ
         return [
             cᴵ,
             rᴵ-cᴵ,
             cᴵ-g,
             rᴵ
         ]
    end
    
    function G_LUI(X)
        return vcat(G_U(X), G_L(X), G_I(X))
    end
    function F_dummy(X)
    end
    function G_dummy(X)
    end

    # four variables
    RSIProblem = MultiLevelProblem(dim)
    n = 10
    if arch == "(LUI)DP"
        RSIProblem.addLevel!(Problem((F_L, F_U, F_I), G_LUI, vcat(rᴸcᴸ_dim, rᵁcᵁ_dim, rᴵ_dim), 5n))
        RSIProblem.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim, Aᴵ_dim), 5n))
        RSIProblem.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p, Pᴵ_dim_p), 20n))
        RSIProblem.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    end

    function RSIProblem.visualize(X; kwargs...)
        rᵁcᵁ = X[rᵁcᵁ_dim]
        rᴸcᴸ = X[rᴸcᴸ_dim]
        rᴵcᴵ = [X[rᴵ_dim][1], X[rᴵ_dim][1] * (1-γ)]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z)
        Pᴾ = reshape(OD - Pᵁ - Pᴸ - Pᴵ, Z, Z)
        println("rᵁcᵁ = ", rᵁcᵁ)
        println("rᴸcᴸ = ", rᴸcᴸ)
        println("rᴵcᴵ = ", rᴵcᴵ)
        println("Aᵁ = ", Aᵁ)
        println("Aᴸ = ", Aᴸ)
        println("Aᴵ = ", Aᴵ)
        println("Pᵁ = ", Pᵁ)
        println("Pᴸ = ", Pᴸ)
        println("Pᴵ = ", Pᴵ)
        println("Pᴾ = ", Pᴾ)
    end

    function find_feasible_point()
        rᵁcᵁ = [2, 1]
        rᴸcᴸ = [1.5, 1]
        rᴵ = [2]
        Aᵁ = Πᴰ ./ 3
        Aᴸ = Πᴰ ./ 3
        Aᴵ = Πᴰ ./ 3
        Pᵁ = OD ./ 4
        Pᴸ = OD ./ 4
        Pᴵ = OD ./ 4
        return vcat(rᵁcᵁ, rᴸcᴸ, rᴵ, Aᵁ, Aᴸ, Aᴵ,
                vec(Pᵁ), vec(Pᴸ), vec(Pᴵ))
    end

    RSIProblem.x_s = find_feasible_point()
    # RSIProblem.x_s = [0.921729805025922, 26.205583958240894, 1.2143599111225507, 4.76329432387886, 25.011495450687793, 1.6877590743896507, 0.4008928153410466, 2.972331763578343, 0.6257519074322478, 0.034245047698801545, 0.2991469250276102, 4.799831936264837, 0.0, 0.10800499029668909, 1.2384515946854138, 0.5273070315738075, 0.0, 1.4099191544988732, 0.7537198708446567, 11.302535013863094, 0.0, 0.0, 1.0395392623228243, 2.4575180050456713, 0.07643642149219798, 0.0, 15.700322251176022, 0.1071785821558191, 0.21377907868414503, 0.0]
    
    X⁰ = [35.58926493624755, 1.0848111260072895, 31.38475653670711, 1.3104034545675232, 0.0026456016159885944, 0.016070883455856045, 0.509392564258502, 0.08884506064769465, 0.3189407829169328, 0.027011096758198197, 0.0, 0.005807456818747853, 0.0014322482124298297, 0.006046606507648412, 0.0, 0.09231895887226814, 0.0001841779184628572, 0.000374818490782235, 0.0, 0.0, 0.02892105112251978, 0.0031552626316130143, 0.016155689339485788, 0.0, 0.0026202236382458386, 0.006194693133560761, 0.03555095214411734, 0.0]
    base_D = -F_D0(X⁰)
    base_P = -F_P0(X⁰)
    
    RSIProblem.MAX_ITER = 150
    RSIProblem.alpha = fill(1., dim)
    RSIProblem.alpha[Aᵁ_dim] .= 0.1
    RSIProblem.alpha[Aᴸ_dim] .= 0.1
    RSIProblem.alpha[Aᴵ_dim] .= 0.1
    RSIProblem.alpha[Pᵁ_dim] .= 0.01
    RSIProblem.alpha[Pᴸ_dim] .= 0.01
    RSIProblem.alpha[Pᴵ_dim] .= 0.01
    export RSIProblem
end