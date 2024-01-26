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
    OD .*= 100      # rescaled so that it works better with learning rate.
    Πᴰ .*= 10

    @assert (sum(Πᴼ) == 1 && sum(Πᵀ) == 1 && sum(Πᴰ) ≈ 10 && sum(OD) ≈ 100)

    # Dimensions
    fᵁrᵁcᵁ_dim = 1:3
    fᴸrᴸcᴸ_dim = 4:6
    fᴵrᴵcᴵ_dim = 7:9
    Aᵁ_dim = 10:(10+Z-1)
    Aᴸ_dim = Aᵁ_dim .+ Z
    Aᴵ_dim = Aᴸ_dim .+ Z
    Pᵁ_dim = collect(Aᴵ_dim[end]+1 : Aᴵ_dim[end]+ Z^2)
    Pᴸ_dim = collect(Pᵁ_dim .+ Z^2)
    Pᴵ_dim = collect(Pᴸ_dim .+ Z^2)

    println("Dimensions:")
    println("fᵁrᵁcᵁ_dim = ", fᵁrᵁcᵁ_dim, " (", length(fᵁrᵁcᵁ_dim), ")")
    println("fᴸrᴸcᴸ_dim = ", fᴸrᴸcᴸ_dim, " (", length(fᴸrᴸcᴸ_dim), ")")
    println("fᴵrᴵcᴵ_dim = ", fᴵrᴵcᴵ_dim, " (", length(fᴵrᴵcᴵ_dim), ")")
    println("Aᵁ_dim = ", Aᵁ_dim, " (", length(Aᵁ_dim), ")")
    println("Aᴸ_dim = ", Aᴸ_dim, " (", length(Aᴸ_dim), ")")
    println("Aᴵ_dim = ", Aᴵ_dim, " (", length(Aᴵ_dim), ")")
    println("Pᵁ_dim = ", Pᵁ_dim, " (", length(Pᵁ_dim), ")")
    println("Pᴸ_dim = ", Pᴸ_dim, " (", length(Pᴸ_dim), ")")
    println("Pᴵ_dim = ", Pᴵ_dim, " (", length(Pᴵ_dim), ")")
    println()
    
    dim = Pᴵ_dim[end]   # Total dimensions

    if Z != 3
        println("Don't forget to ammend the snippet below")
    end
    # Diagonal elements in demand matrices are always 0 (and unchanged)
    Pᵁ_dim_p = setdiff(Pᵁ_dim, Pᵁ_dim[1:Z+1:end])
    Pᴸ_dim_p = setdiff(Pᴸ_dim, Pᴸ_dim[1:Z+1:end])
    Pᴵ_dim_p = setdiff(Pᴵ_dim, Pᴵ_dim[1:Z+1:end])
    
    Aᴾ = [10, 10, 10]  # availability of public transit
    λ = 0.5
    g = 1               # gas cost
    fᴾ = 15             # fixed initial cost for public vehicles
    rᴾ = 0.5            # rates for public vehicles

    function F_P(X)         # Passenger's objective
        ϵ = 1e-6
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z)
        Pᴾ = reshape(OD - Pᵁ - Pᴸ - Pᴵ, Z, Z)
        fᵁ, rᵁ, _ = X[fᵁrᵁcᵁ_dim]
        fᴸ, rᴸ, _ = X[fᴸrᴸcᴸ_dim]
        fᴵ, rᴵ, _ = X[fᴵrᴵcᴵ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]
        Pᵁᵢ = sum(Pᵁ, dims=2)
        Pᴸᵢ = sum(Pᴸ, dims=2)
        Pᴵᵢ = sum(Pᴵ, dims=2)
        Pᴾᵢ = sum(Pᴾ, dims=2)
        acc = 0.
        for i in 1:Z
            for j in 1:Z
                acc += Pᵁ[i,j] * (fᵁ + rᵁ * C[i, j] + λ * Pᵁᵢ[i]/(Aᵁ[i]+ϵ)) +
                       Pᴸ[i,j] * (fᴸ + rᴸ * C[i, j] + λ * Pᴸᵢ[i]/(Aᴸ[i]+ϵ)) +
                       Pᴵ[i,j] * (fᴵ + rᴵ * C[i, j] + λ * Pᴵᵢ[i]/(Aᴵ[i]+ϵ)) +
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
        return vcat(vec(Pᵁ), vec(Pᴸ), vec(Pᴵ), vec(Pᴾ))
    end

    function F_D(X)         # Driver's objective
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z)
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z)
        # Pᴾ = reshape(OD - Pᵁ - Pᴸ - Pᴵ, Z, Z)
        _, _, cᵁ = X[fᵁrᵁcᵁ_dim]
        _, _, cᴸ = X[fᴸrᴸcᴸ_dim]
        _, _, cᴵ = X[fᴵrᴵcᴵ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]
        acc = 0
        for i in 1:Z
            acc += Aᵁ[i] * (cᵁ - g) * sum(Pᵁ[i, j]*C[i, j] for j in 1:Z) + 
                    Aᴸ[i] * (cᴸ - g) * sum(Pᴸ[i,j]*C[i, j] for j in 1:Z) + 
                    Aᴵ[i] * (cᴵ - g) * sum(Pᴵ[i,j]*C[i, j] for j in 1:Z)
        end
        return -acc     # maximization
    end

    function G_D(X)
        _, _, cᵁ = X[fᵁrᵁcᵁ_dim]
        _, _, cᴸ = X[fᴸrᴸcᴸ_dim]
        _, _, cᴵ = X[fᴵrᴵcᴵ_dim]
        Aᵁ = X[Aᵁ_dim]
        Aᴸ = X[Aᴸ_dim]
        Aᴵ = X[Aᴵ_dim]
        return vcat(Πᴰ - Aᵁ - Aᴸ - Aᴵ, Aᵁ, Aᴸ, Aᴵ, [cᵁ-g, cᴸ-g, cᴵ-g])
    end

    function F_U(X)         # Uber
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z)
        fᵁ, rᵁ, cᵁ = X[fᵁrᵁcᵁ_dim]
        acc = 0.
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

    function G_LUI(X)
        return vcat(G_U(X), G_L(X), G_I(X))
    end

    function F_I(X)
       util_D = -F_D(X) # utility for driver 
       util_P = -F_P(X) # utility for passenger
       return -(util_D - base_D)*(util_P - base_P)  # max
    end

    function G_I(X)
        fᴵ, rᴵ, cᴵ = X[fᴵrᴵcᴵ_dim]
        return [
            cᴵ,
            cᴵ-g,
            fᴵ,
            rᴵ
        ]
    end
    
    function F_dummy(X)
    end
    function G_dummy(X)
    end

    # four variables
    RSProblemI = MultiLevelProblem(dim)
    n = 10
    if arch == "(LUI)DP"
        RSProblemI.addLevel!(Problem((F_L, F_U, F_I), G_LUI, vcat(fᴸrᴸcᴸ_dim, fᵁrᵁcᵁ_dim, fᴵrᴵcᴵ_dim), 5n))
        RSProblemI.addLevel!(Problem(F_D, G_D, vcat(Aᵁ_dim, Aᴸ_dim, Aᴵ_dim), 5n))
        RSProblemI.addLevel!(Problem(F_P, G_P, vcat(Pᵁ_dim_p, Pᴸ_dim_p, Pᴵ_dim_p), 20n))
        RSProblemI.addLevel!(Problem(F_dummy, G_dummy, [], 0))   # dummy won't have any effect
    end

    function RSProblemI.visualize(X; kwargs...)
        fᵁrᵁcᵁ = X[fᵁrᵁcᵁ_dim]
        fᴸrᴸcᴸ = X[fᴸrᴸcᴸ_dim]
        fᴵrᴵcᴵ = X[fᴵrᴵcᴵ_dim]
        Aᵁ = X[Aᵁ_dim] ./ 10
        Aᴸ = X[Aᴸ_dim] ./ 10
        Aᴵ = X[Aᴵ_dim] ./ 10
        Pᵁ = reshape(X[Pᵁ_dim], Z, Z) ./ 100
        Pᴸ = reshape(X[Pᴸ_dim], Z, Z) ./ 100
        Pᴵ = reshape(X[Pᴵ_dim], Z, Z) ./ 100
        Pᴾ = Pᴾ = reshape(OD - Pᵁ - Pᴸ - Pᴵ, Z, Z) ./ 100
        println("fᵁrᵁcᵁ = ", fᵁrᵁcᵁ)
        println("fᴸrᴸcᴸ = ", fᴸrᴸcᴸ)
        println("fᴸrᴸcᴸ = ", fᴵrᴵcᴵ)
        println("Aᵁ = ", Aᵁ)
        println("Aᴸ = ", Aᴸ)
        println("Aᴵ = ", Aᴵ)
        println("Pᵁ = ", Pᵁ)
        println("Pᴸ = ", Pᴸ)
        println("Pᴵ = ", Pᴵ)
        println("Pᴾ = ", Pᴾ)
    end

    function find_feasible_point()
        fᵁrᵁcᵁ = [5, 2, 1]
        fᴸrᴸcᴸ = [4, 1.5, 1]
        fᴵrᴵcᴵ = [4, 1, 1]
        Aᵁ = Πᴰ ./ 3
        Aᴸ = Πᴰ ./ 3
        Aᴵ = Πᴰ ./ 3
        Pᵁ = OD ./ 4
        Pᴸ = OD ./ 4
        Pᴵ = OD ./ 4
        return vcat(fᵁrᵁcᵁ, fᴸrᴸcᴸ, fᴵrᴵcᴵ, Aᵁ, Aᴸ, Aᴵ,
                vec(Pᵁ), vec(Pᴸ), vec(Pᴵ))
    end

    X⁰ = [16.294139789326557, 40.45474298175821, 1.4463177264781892, 7.541067156629865, 31.81584493605669, 1.0208981453894346, 0.04685972790665238, 0.16040302203818446, 5.29925933411386, 0.33650012379509864, 2.9476355802374825, 3.256501692894265e-5, 0.0, 0.6267947792749295, 0.49071632940784804, 0.3226401661977072, 0.0, 13.684000325509942, 0.049669551183246785, 0.9823351452150004, 0.0, 0.0, 2.282943552025342, 0.019618851101434576, 0.4312772555534341, 0.0, 0.0045658670730714945, 0.7432887431486215, 6.337985951287279, 0.0]
    base_D = -F_D0(X⁰)
    base_P = -F_P0(X⁰)
    println("Baseline for Driver, Passenger", base_D, base_P)
    RSProblemI.x_s = find_feasible_point()
    # RSProblem.x_s = [0.921729805025922, 26.205583958240894, 1.2143599111225507, 4.76329432387886, 25.011495450687793, 1.6877590743896507, 0.4008928153410466, 2.972331763578343, 0.6257519074322478, 0.034245047698801545, 0.2991469250276102, 4.799831936264837, 0.0, 0.10800499029668909, 1.2384515946854138, 0.5273070315738075, 0.0, 1.4099191544988732, 0.7537198708446567, 11.302535013863094, 0.0, 0.0, 1.0395392623228243, 2.4575180050456713, 0.07643642149219798, 0.0, 15.700322251176022, 0.1071785821558191, 0.21377907868414503, 0.0]
    RSProblemI.MAX_ITER = 150
    RSProblemI.alpha = 1
    export RSProblemI
end