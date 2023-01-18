module Bandsolvers

using LinearAlgebra: eigen, schur, diagm, diagind, eigvals, ⋅, mul!, Diagonal, Symmetric, Hermitian

"A type for storing the Wannier functions."
mutable struct Wanniers
    minlevel::Int # number of the first energy level of ℎ to use for constructing wanniers (this is used in the unperturbed case)
    targetlevels::Vector{Int} # numbers of quasienergy levels to use for constructing wanniers (this is used in the Floquet case)
    n_lo::Vector{Int} # number of levels in the lower subband at each phase; in non-periodic case this depends on the edge state branch position 
    E::Matrix{Float64} # `E[j, i]` = mean energy of `j`th wannier at `i`th phase
    pos::Matrix{Float64} # `pos[j, i]` = position eigenvalue of `j`th wannier at `i`th phase
    d::Array{ComplexF64, 3} # `d[:, :, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details
end

"Default-construct an empty `Wanniers` object."
Wanniers() = Wanniers(0, Int[], Int[], Float64[;;], Float64[;;], ComplexF64[;;;])

"""
A type representing the unperturbed Hamiltonian
    ℎ = 𝑝²/2𝑀 + 𝑔ₗcos²(2𝑥) + 𝑉ₗcos²(𝑥 + 𝜑ₓ).
"""
mutable struct UnperturbedHamiltonian
    N::Int  # number of lattice cells
    M::Float64
    l::Int
    gₗ::Float64
    Vₗ::Float64
    isperiodic::Bool
    φₓ::Vector{Float64}
    maxlevel::Int   # highest level number to consider
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, 1 ≤ i ≤ maxlevel, 1 ≤ j ≤ length(φₓ)
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::Wanniers
end

"""
Construct an `UnperturbedHamiltonian` object. `maxband` is the highest energy band number to consider.
Each band is assumed to consist of 2 subbands containing `2n_cells` levels in total.
"""
function UnperturbedHamiltonian(n_cells::Integer; M::Real, gₗ::Real, Vₗ::Real, maxband::Integer, isperiodic::Bool, φₓ::AbstractVector{<:Real},
                                l::Union{Nothing, Integer}=nothing)
    # convert max band number to level number
    if isperiodic
        maxlevel = maxband * 2n_cells
    else
        maxlevel = (maxband-1) ÷ 2 * 4n_cells + (isodd(maxband) ? 2n_cells-1 : 4n_cells)
    end

    E = Matrix{Float64}(undef, maxlevel, length(φₓ))
    c = Array{ComplexF64,3}(undef, 2maxlevel+1, maxlevel, length(φₓ))

    UnperturbedHamiltonian(Int(n_cells), Float64(M), (l === nothing ? 1 : l), Float64(gₗ), Float64(Vₗ), isperiodic,
                           collect(Float64, φₓ), maxlevel, E, c, Wanniers())
end

"Diagonalise the unperturbed Hamiltonian `uh` at each phase."
function diagonalise!(uh::UnperturbedHamiltonian)
    (;N, M, gₗ, Vₗ, maxlevel) = uh
    sortby = M > 0 ? (+) : (-) # eigenvalue sorting; for 𝑀 < 0 we use descending sorting
    if uh.isperiodic
        h = diagm(0 => ComplexF64[(2j/N)^2 / 2M + (gₗ + Vₗ)/2 for j = -maxlevel:maxlevel])
        h[diagind(h, 2N)] .= gₗ/4
        for (i, φ) in enumerate(uh.φₓ)
            h[diagind(h, N)] .= Vₗ/4 * cis(-2φ)
            f = eigen(Hermitian(h); sortby)
            uh.E[:, i] = f.values[1:maxlevel]
            uh.c[:, :, i] = f.vectors[:, 1:maxlevel]
        end
    else
        X(j′, j) = 16N*j*j′ / (π*((j-j′)^2-(2N)^2)*((j+j′)^2-(2N)^2))
        n_j = 2maxlevel + 1
        h = zeros(n_j, n_j)
        for (i, φ) in enumerate(uh.φₓ)
            for j in 1:n_j
                for j′ in 1:j
                    val = 0.0
                    if isodd(j′ + j)
                        val += Vₗ/2 * X(j′, j) * sin(2φ)
                    else
                        # check diagonals "\"
                        if j′ == j
                            val += (gₗ + Vₗ)/2 + (j / N)^2 / 2M
                        elseif j′ == j - 2N
                            val += Vₗ * cos(2φ) / 4
                        elseif j′ == j - 4N
                            val += gₗ/4
                        end
                        # check anti-diagonals "/"
                        if j′ == -j + 2N
                            val += -Vₗ * cos(2φ) / 4
                        elseif j′ == -j + 4N
                            val += -gₗ/4
                        end
                    end
                    h[j′, j] = val
                end
            end
            f = eigen(Symmetric(h); sortby)
            uh.E[:, i] = f.values[1:maxlevel]
            uh.c[:, :, i] = f.vectors[:, 1:maxlevel]
        end
    end
end

"""
Calculate Wannier vectors for the unperturbed Hamiltonian using the energy eigenstates in the band number `targetband`.
`mixsubbands` indicates whether the two subbands have to be mixed or treated separately.
"""
function compute_wanniers!(uh::UnperturbedHamiltonian; targetband::Integer, mixsubbands::Bool)
    N = uh.N

    minlevel = (targetband-1) * 2N + 1

    if uh.isperiodic
        E = Matrix{Float64}(undef, 2N, length(uh.φₓ))
        pos = Matrix{Float64}(undef, 2N, length(uh.φₓ))
        
        if mixsubbands
            d = Array{ComplexF64, 3}(undef, 2N, 2N, length(uh.φₓ))
            n_lo = fill(2N, length(uh.φₓ)) # treat all Wanniers as belonging to the lower subband when plotting (see `make_wannierfunctions`)
            uh.w = Wanniers(minlevel, Int[], n_lo, E, pos, d)

            X = Matrix{ComplexF64}(undef, 2N, 2N) # position operator
            for i in eachindex(uh.φₓ)
                for n in 1:2N
                    for n′ in 1:2N
                        @views X[n′, n] = uh.c[2:end, minlevel+n′-1, i] ⋅ uh.c[1:end-1, minlevel+n-1, i]
                    end
                end
                _, d[:, :, i], pos_complex = schur(X)
                pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
                sp = sortperm(pos_real)                 # sort the eigenvalues
                pos[:, i] = pos_real[sp]
                @views Base.permutecols!!(d[:, :, i], sp)    # sort the eigenvectors in the same way
                E[:, i] = transpose(uh.E[range(minlevel, length=2N), i]) * abs2.(d[:, :, i])
            end
        else
            # `d` fill format: `d[1:N, 1:N, i]` = eigenvectors of the lower subband,
            #                  `d[1:N, N+1:2N, i]` = eigenvectors of the higher subband
            d = Array{ComplexF64, 3}(undef, N, 2N, length(uh.φₓ))
            uh.w = Wanniers(minlevel, Int[], fill(N, length(uh.φₓ)), E, pos, d)

            X = Matrix{ComplexF64}(undef, N, N) # position operator
            for i in eachindex(uh.φₓ)
                for o in (0, N)
                    window = 1+o:N+o
                    for n in 1:N
                        for n′ in 1:N
                            @views X[n′, n] = uh.c[2:end, minlevel+o+n′-1, i] ⋅ uh.c[1:end-1, minlevel+o+n-1, i]
                        end
                    end
                    _, d[:, window, i], pos_complex = schur(X)
                    pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
                    sp = sortperm(pos_real)                 # sort the eigenvalues
                    pos[window, i] = pos_real[sp]
                    @views Base.permutecols!!(d[:, window, i], sp)    # sort the eigenvectors in the same way
                    E[window, i] = transpose(uh.E[range(minlevel+o, length=N), i]) * abs2.(d[:, window, i])
                end
            end
        end
    else # if !uh.isperiodic
        n_w = isodd(targetband) ? 2N-1 : 2N+1 # total number of wanniers to construct; this is the number of levels in the target band
        E = Matrix{Float64}(undef, n_w, length(uh.φₓ))
        pos = Matrix{Float64}(undef, n_w, length(uh.φₓ))
        
        n_j = size(uh.c, 1)
        
        if mixsubbands
            d = Array{ComplexF64, 3}(undef, n_w, n_w, length(uh.φₓ))
            uh.w = Wanniers(minlevel, Int[], fill(n_w, length(uh.φₓ)), E, pos, d)

            X = Matrix{ComplexF64}(undef, n_w, n_w) # position operator
            for i in eachindex(uh.φₓ)
                for n in 1:n_w
                    for n′ in n:n_w
                        X[n′, n] = X[n, n′] = (n == n′ ? N*π/2 : 0.0) - 8N/π*sum(uh.c[j, minlevel+n-1, i] * sum(uh.c[j′, minlevel+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                                 for j′ = (iseven(j) ? 1 : 2):2:n_j) for j = 1:n_j)
                    end
                end
                uh.w.pos[:, i], uh.w.d[:, :, i] = eigen(Hermitian(X))
                uh.w.E[:, i] = transpose(uh.E[range(minlevel, length=n_w), i]) * abs2.(uh.w.d[:, :, i])
            end
        else
            # `d` fill format: `d[1:n_lo[i], 1:n_lo[i], i]` = eigenvectors of the lower subband,
            #                  `d[1:n_w-n_lo[i], n_lo[i]+1:n_w, i]` = eigenvectors of the higher subband
            d = Array{ComplexF64, 3}(undef, n_w÷2+1, n_w, length(uh.φₓ))
            uh.w = Wanniers(minlevel, Int[], Vector{Int}(undef, length(uh.φₓ)), E, pos, d)
            
            X = ComplexF64[;;] # position operator
            X_less = zeros(n_w÷2, n_w÷2) # position operator for a subband which does not contain the edge state branch
            X_more = zeros(n_w÷2+1, n_w÷2+1) # position operator for a subband which contains the edge state branch
            for i in eachindex(uh.φₓ)
                up = uh.E[minlevel+n_w÷2, i] > (uh.E[minlevel+n_w÷2-1, i] + uh.E[minlevel+n_w÷2+1, i])/2 # true if the edge state branch is above the mean value
                
                # Lower band
                X = up ? X_less : X_more # bind the position operator to a matrix of the appropriate size
                n_lo = n_w÷2 + !up # number of levels in the lower subband
                uh.w.n_lo[i] = n_lo
                for n in 1:n_lo
                    for n′ in n:n_lo
                        X[n′, n] = X[n, n′] = (n == n′ ? N*π/2 : 0.0) - 8N/π*sum(uh.c[j, minlevel+n-1, i] * sum(uh.c[j′, minlevel+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                                 for j′ = (iseven(j) ? 1 : 2):2:n_j) for j = 1:n_j)
                    end
                end
                uh.w.pos[1:n_lo, i], uh.w.d[1:n_lo, 1:n_lo, i] = eigen(Hermitian(X))
                uh.w.E[1:n_lo, i] = transpose(uh.E[range(minlevel, length=n_lo), i]) * abs2.(uh.w.d[1:n_lo, 1:n_lo, i])

                # Higher band
                X = up ? X_more : X_less
                n_hi = n_w - n_lo
                for n in 1:n_hi
                    for n′ in n:n_hi
                        X[n′, n] = X[n, n′] = (n == n′ ? N*π/2 : 0.0) - 8N/π*sum(uh.c[j, minlevel+n_lo+n-1, i] * sum(uh.c[j′, minlevel+n_lo+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                                 for j′ = (iseven(j) ? 1 : 2):2:n_j) for j = 1:n_j)
                    end
                end
                uh.w.pos[n_lo+1:n_w, i], uh.w.d[1:n_hi, n_lo+1:n_w, i] = eigen(Hermitian(X))
                uh.w.E[n_lo+1:n_w, i] = transpose(uh.E[range(minlevel+n_lo, length=n_hi), i]) * abs2.(uh.w.d[1:n_hi, n_lo+1:n_w, i])
            end
        end
    end
end

"""
Construct energy eigenfunctions at coordinates in `x` for each eigenstate number in `whichstates` at each phase number in `whichphases`.
Return `ψ`, where `ψ[:, j, i]` = `j`th eigenfunction at `whichphases[i]`th phase.
"""
function make_eigenfunctions(uh::UnperturbedHamiltonian, x::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer}, whichstates::AbstractVector{<:Integer})
    ψ = Array{ComplexF64,3}(undef, length(x), length(whichstates), length(whichphases))
    make_state = uh.isperiodic ? make_exp_state : make_sin_state
    for (i, iφ) in enumerate(whichphases)
        for (j, js) in enumerate(whichstates)
            @views ψ[:, j, i] = make_state(x, uh.c[:, js, iφ]; N=uh.N)
        end
    end
    return ψ
end

"""
Construct Wannier functions at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `uh` are constructed.
In the process, energy eigenfunctions are also constructed.
Return `ψ, w`, where `ψ[:, j, i]` = `j`th eigenfunction at `whichphases[i]`th phase, and `w[:, j, i]` = `j`th Wannier function at `i`th phase.
"""
function make_wannierfunctions(uh::UnperturbedHamiltonian, x::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
    n_w = size(uh.w.E, 1)
    w = Array{ComplexF64, 3}(undef, length(x), n_w, length(whichphases))
    ψ = make_eigenfunctions(uh, x, whichphases, range(uh.w.minlevel, length=n_w))
    for (i, iφ) in enumerate(whichphases)
        window = 1:uh.w.n_lo[i]
        w[:, window, i] = ψ[:, window, i] * uh.w.d[window, window, iφ]
        window = uh.w.n_lo[i]+1:n_w 
        window2 = 1:n_w-uh.w.n_lo[i]
        w[:, window, i] = ψ[:, window, i] * uh.w.d[window2, window, iφ]
    end
    return ψ, w
end

"Construct the coordinate-space wavefunction 𝜓(𝑥) = ∑ⱼ𝑐ⱼexp(2i𝑗𝑥/𝑁) / √(𝑁π)"
function make_exp_state(x::AbstractVector{<:Real}, c::AbstractVector{<:Number}; N)
    ψ = zeros(eltype(c), length(x))
    n_j = (length(c) - 1) ÷ 2
    for j in -n_j:n_j
        @. ψ += c[j+n_j+1] * cis(2j/N * x)
    end
    return ψ ./ sqrt(N*π)
end

"Construct the coordinate-space wavefunction 𝜓(𝑥) = ∑ⱼ𝑐ⱼsin(𝑗𝑥/𝑁) / √(𝑁π/2)"
function make_sin_state(x::AbstractVector{<:Real}, c::AbstractVector{<:Number}; N)
    ψ = zeros(eltype(c), length(x))
    for (j, c) in enumerate(c)
        @. ψ += c * sin(j/N * x)
    end
    return ψ ./ sqrt(N*π/2)
end

"""
A type representing a tight-binding Hamiltonian.
"""
mutable struct TBHamiltonian
    N::Int
    H::Array{ComplexF64, 3} # Hamiltonian matrix
    φₓ::Vector{Float64}
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, `1 ≤ i ≤ 2N`, `1 ≤ j ≤ length(φₓ)`
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::Wanniers 
end

"Construct a `TBHamiltonian` object. `uh` must contain calculated periodic Wanniers."
function TBHamiltonian(uh::UnperturbedHamiltonian)
    (;N, gₗ, Vₗ, φₓ) = uh
    n_φₓ = length(φₓ)
    n_w = size(uh.w.E, 1) # number of Wanniers
    H = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ) # TB Hamiltonian matrix
    iφ₀ = 1 # phase index at which to take the Wanniers -- any choice is OK
    # Matrix of Wannier basis vectors |𝑤ₐ⟩ = ∑ᵢ 𝑑ᵃᵢ |𝜓ᵢ⟩
    w = uh.c[:, range(uh.w.minlevel, length=n_w), iφ₀] * uh.w.d[:, :, iφ₀]
    
    # Matrix of the unperturbed Hamiltonian
    h = diagm(0 => ComplexF64[(2j/N)^2 / 2uh.M + (gₗ + Vₗ)/2 for j = -uh.maxlevel:uh.maxlevel])
    h[diagind(h, -2N)] .= h[diagind(h, 2N)] .= gₗ/4
    # Compute elements of `H` at each phase
    for (iφ, φ) in enumerate(φₓ)
        h[diagind(h, -N)] .= Vₗ/4 * cis(+2φ)
        h[diagind(h, +N)] .= Vₗ/4 * cis(-2φ)
        H[:, :, iφ] = w' * h * w
    end
    E = Matrix{Float64}(undef, n_w, n_φₓ)
    c = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ)
    TBHamiltonian(N, H, φₓ, E, c, Wanniers())
end

"Diagonalise the TB Hamiltonian `tbh` at each phase."
function diagonalise!(tbh::TBHamiltonian)
    for iφ in eachindex(tbh.φₓ)
        tbh.E[:, iφ], tbh.c[:, :, iφ] = eigen(Hermitian(tbh.H[:, :, iφ]))
    end
end

"Calculate Wannier vectors for each of the two subbands for the TB Hamiltonian `tbh`."
function compute_wanniers!(tbh::TBHamiltonian)
    (;N, φₓ) = tbh
    # `d` fill format: `d[1:N, 1:N, i]` = eigenvectors of the lower subband,
    #                  `d[1:N, N+1:2N, i]` = eigenvectors of the higher subband
    d = Array{ComplexF64, 3}(undef, N, 2N, length(φₓ))
    E = Matrix{Float64}(undef, 2N, length(φₓ))
    pos = Matrix{Float64}(undef, 2N, length(φₓ))
    tbh.w = Wanniers(0, Int[], Int[], E, pos, d)

    X = Diagonal([cis(2π/(N*π) * n*π/2) for n in 0:2N-1]) # position operator in coordinate representation
    for b in 1:2
        levels = N*(b-1)+1:N*b
        for iφ in eachindex(φₓ)
            XE = tbh.c[:, levels, iφ]' * X * tbh.c[:, levels, iφ] # position operator in energy representation
            _, d[:, levels, iφ], pos_complex = schur(XE)
            pos_real = @. mod2pi(angle(pos_complex)) / 2π * N*π # shift angle from [-π, π) to [0, 2π)
            sp = sortperm(pos_real)                        # sort the eigenvalues
            pos[levels, iφ] = pos_real[sp]
            @views Base.permutecols!!(d[:, levels, iφ], sp) # sort the eigenvectors in the same way
            E[levels, iφ] = transpose(tbh.E[levels, iφ]) * abs2.(d[:, levels, iφ])
        end
    end
end

"""
Construct Wannier functions at each phase number in `whichphases`. All Wannier functions contained in `h` are constructed.
Return `w`, where `w[:, j, i]` = `j`th Wannier function at `i`th phase.
"""
function make_wannierfunctions(tbh::TBHamiltonian, whichphases::AbstractVector{<:Integer})
    (;N) = tbh
    n_w = size(tbh.w.E, 1)
    w = Array{ComplexF64, 3}(undef, size(tbh.c, 1), n_w, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for b in 1:2
            levels = N*(b-1)+1:N*b
            w[:, levels, i] = tbh.c[:, range(N*(b-1)+1, length=N), iφ] * tbh.w.d[:, levels, iφ]
        end
    end
    return w
end

"""
A type representing the Floquet Hamiltonian
    ℋ = ℎ - i∂ₜ + λₛsin²(2𝑥)cos(2𝜔𝑡) + λₗcos²(2𝑥)cos(𝜔𝑡 + 𝜑ₜ),
where ℎ is the unperturbed Hamiltonian represented by [`UnperturbedHamiltonian`](@ref), and 𝜑ₜ = 2𝜑ₓ.
"""
mutable struct FloquetHamiltonian
    uh::UnperturbedHamiltonian
    s::Int
    λₛ::Float64
    λₗ::Float64
    ω::Float64
    pumptype::Symbol
    minlevel::Int # lowest level number of ℎ to use when constructing ℋ
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue (Floquet quasienergy) at `j`th phase, `i = 1:maxlevel`
    b::Array{ComplexF64, 3} # `b[:, i, j]` = `i`th eigenvector at `j`th phase, `i = 1:maxlevel; j = 1:2maxlevel+1`
    ν::Vector{Int}  # band map 𝜈(𝑚)
end

"""
Construct a `FloquetHamiltonian` object. `minband` is the first energy band of `uh` to use when constructing the Floquet Hamiltonian matrix.
Type of pumping is controlled via `pumptype`: `:time` for temporal, `:space` for spatial, or anything else for simultaneous space-time pumping.
In the case of time-only pumping, it is assumed that 𝜑ₓ = 0, and hence that `uh.φₓ[1] == 0`.
"""
function FloquetHamiltonian(uh::UnperturbedHamiltonian; s::Integer, λₛ::Real, λₗ::Real, ω::Real, pumptype::Symbol, minband::Integer)
    N = uh.N

    bs1 = 2N - 1 # "bandsize 1" = number of levels in the lowest band of ℎ

    # convert min band number to level number
    if uh.isperiodic
        minlevel = (minband - 1) * 2N + 1
    else
        minlevel = (minband - 1) ÷ 2 * 4N + (isodd(minband) ? 1 : bs1 + 1)
    end
    
    if uh.isperiodic
        ν = [ceil(Int, m/2N) for m in 1:uh.maxlevel] 
    else
        # FIll `ν`: [1 (`bs1` times), 2 (`bs2 = 2N+1` times), 3 (`bs1` times), 4 (`bs2` times), ...]
        ν = Vector{Int}(undef, uh.maxlevel)
        number = 1
        for i in 0:uh.maxlevel÷4N-1
            ν[4N*i+1:4N*i+bs1] .= number
            number += 1
            ν[4N*i+bs1+1:4N*i+4N] .= number
            number += 1
        end
        ν[uh.maxlevel - uh.maxlevel%4N + 1:end] .= number
    end
    
    n_levels = uh.maxlevel - minlevel + 1
    E = Matrix{Float64}(undef, n_levels, length(uh.φₓ))
    b = Array{ComplexF64,3}(undef, n_levels, n_levels, length(uh.φₓ))
    
    FloquetHamiltonian(uh, Int(s), Float64(λₛ), Float64(λₗ), Float64(ω), pumptype, Int(minlevel), E, b, ν)
end

"Diagonalise the Floquet Hamiltonian `fh` at each phase."
function diagonalise!(fh::FloquetHamiltonian)
    (;N, φₓ, E, c) = fh.uh
    n_j = size(c, 1)
    (;s, ω, λₛ, λₗ, pumptype, ν) = fh

    n_levels = fh.uh.maxlevel - fh.minlevel + 1 # number of levels of spatial Hamiltonian to use for constructing ℋ

    H = zeros(ComplexF64, n_levels, n_levels) # ℋ matrix, will only fill the lower triangle

    ∑cc(m′, m, i) = if fh.uh.isperiodic
        @views c[1+2N:end, m′, i] ⋅ c[1:end-2N, m, i] + c[1:end-2N, m′, i] ⋅ c[1+2N:end, m, i]
    else
        sum( (-c[-j+4N, m′, i] + c[j+4N, m′, i]) * c[j, m, i] for j = 1:4N-1 ) +
                               + c[4N+4N, m′, i] * c[4N, m, i] + # iteration `j = 4N`
        sum( (  c[j-4N, m′, i] + c[j+4N, m′, i]) * c[j, m, i] for j = 4N+1:n_j-4N ) + 
        sum( (  c[j-4N, m′, i]                 ) * c[j, m, i] for j = n_j-4N+1:n_j )
    end

    if fh.uh.isperiodic
        for (i, φ) in enumerate(φₓ)
            # `m` and `m′` number the levels of ℎ
            # `e` and `e′` number the elements of `H`
            for m in fh.minlevel:fh.uh.maxlevel
                e = m - fh.minlevel + 1

                # for time-only pumping, always take the eigenenergies at the first phase, which is asssumed to correspond to 𝜑ₓ = 0
                p = (pumptype == :time ? 1 : i)
                H[e, e] = E[m, p] - ν[m]*ω/s

                # place the elements of the long lattice
                for g in 1:2N
                    m′ = 2N*(s + ν[m] - 1) + g
                    e′ = m′ - fh.minlevel + 1
                    e′ > n_levels && break
                    if pumptype != :time || i == 1 # if pumping is time-only, this must be calculated only once, at `i` = 1
                        # if pumping is space-time, then also multiply by cis(-𝜑ₜ). `φ` runs over 𝜑ₓ, and we assume the pumping protocol 𝜑ₜ = 2𝜑ₓ
                        H[e′, e] = (pumptype == :space ? λₗ/8 * ∑cc(m′, m, i) : λₗ/8 * ∑cc(m′, m, i) * cis(-2φ))
                    elseif pumptype == :time 
                        H[e′, e] *= cis(-2(φₓ[i]-φₓ[i-1]))
                    end
                end
                
                # place the elements of the short lattice
                if pumptype != :time || i == 1 # if pumping is time-only, this must be calculated only once, at `i` = 1
                    for g in 1:2N
                        m′ = 2N*(2s + ν[m] - 1) + g
                        e′ = m′ - fh.minlevel + 1
                        e′ > n_levels && break
                        H[e′, e] = -λₛ/8 * ∑cc(m′, m, i)
                    end
                end
            end
            fh.E[:, i], fh.b[:, :, i] = eigen(Hermitian(H, :L), sortby=-)
        end
    else
        for (i, φ) in enumerate(φₓ)
            bs1 = 2N - 1
            bs2 = 2N + 1
            pattern = [fill(bs1, bs1); fill(bs2, bs2)]
            G = repeat(pattern, fh.uh.maxlevel÷4N) # a pattern which e.g. for `N == 2` looks like [3, 3, 3, 5, 5, 5, 5, 3, 3, 3, 5, 5, 5, 5, ...]
            fh.uh.maxlevel % 4N != 0 && append!(G, fill(bs1, bs1))

            # `m` and `m′` number the levels of ℎ
            # `e` and `e′` number the elements of `H`
            for m in fh.minlevel:fh.uh.maxlevel
                e = m - fh.minlevel + 1

                # for time-only pumping, always take the eigenenergies at the first phase, corresponding to 𝜑ₓ = 0
                p = (pumptype == :time ? 1 : i)
                H[e, e] = E[m, p] - fh.ν[m]*ω/s

                # place the elements of the long lattice
                for g in 1:G[m]
                    # skip `s` groups of `4N`, then some more groups depending on `m`, then skip `G[fh.minlevel]` cells
                    e′ = 4N*(s÷2) + 4N*((fh.ν[m]-1)÷2) + iseven(ν[m])*G[fh.minlevel] + g
                    e′ > n_levels && break
                    m′ = e′ + fh.minlevel - 1 
                    if pumptype != :time || i == 1 # if pumping is time-only, this must be calculated only once, at `i` = 1
                        # if pumping is space-time, then also multiply by cis(-𝜑ₜ). `φ` runs over 𝜑ₓ, and we assume the pumping protocol 𝜑ₜ = 2𝜑ₓ
                        H[e′, e] = (pumptype == :space ? λₗ/8 * ∑cc(m′, m, i) : λₗ/8 * ∑cc(m′, m, i) * cis(-2φ))
                    elseif pumptype == :time 
                        H[e′, e] *= cis(-2(φₓ[i]-φₓ[i-1]))
                    end
                end
                
                # place the elements of the short lattice
                if pumptype != :time || i == 1 # if pumping is time-only, this must be calculated only once, at `i` = 1
                    for g in 1:2N
                        e′ = 4N*s + 4N*((fh.ν[m]-1)÷2) + iseven(ν[m])*G[fh.minlevel] + g
                        e′ > n_levels && break
                        m′ = e′ + fh.minlevel - 1
                        H[e′, e] = -λₛ/8 * ∑cc(m′, m, i)
                    end
                end
            end
            fh.E[:, i], fh.b[:, :, i] =  eigen(Hermitian(H, :L), sortby=-)
        end
    end
end

"""
Construct Floquet modes at coordinates in `x` and time moments in `Ωt` for each state number in `whichstates` at each phase number in `whichphases`.
Return `u`, where `u[ix, it, j, i]` = `j`th wavefunction at `whichphases[i]`th phase at `ix`th coordinate at `it`th time moment.
"""
function make_eigenfunctions(fh::FloquetHamiltonian, x::AbstractVector{<:Real}, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer},
                             whichstates::AbstractVector{<:Integer})
    u = Array{ComplexF64,4}(undef, length(x), length(Ωt), length(whichstates), length(whichphases))
    n_levels = size(fh.E, 1) # number of Floquet levels; equivalently, number of levels of ℎ used for constructing ℋ
    # Eigenfunctions of ℎ, which are mixed during construction of `u`. For time-only pumping use only eigenstates at the first phase, corresponding to 𝜑ₓ = 0
    ψ = make_eigenfunctions(fh.uh, x, (fh.pumptype == :time ? [1] : whichphases), range(fh.minlevel, length=n_levels))
    for (i, iφ) in enumerate(whichphases)
        p = (fh.pumptype == :time ? 1 : i)
        for (j, js) in enumerate(whichstates)
            for (it, t) in enumerate(Ωt)
                u[:, it, j, i] = sum(cis(-fh.ν[m]*t) * ψ[:, m, p] * fh.b[m, js, iφ] for m in 1:n_levels)
            end
        end
    end
    return u
end

"""
Permute Floquet quasienergy levels contained in `fh.E` so that they are stored in the same order as the eigenenergies of ℎ stored in `fh.uh.E`.
Repeat this for every phase (i.e. column of `fh.E`).
To perfrorm the sorting, first calculate `fh.uh.E - fh.ν[m]`, which is the diagonal of ℋ. If there is no perturbation, then these
are the Floquet quasienergies. Then, sort them in descending order (as if we diagonalised the Hamiltonian) and find the permutation
that would undo this sorting. This permutation is applied to a copy of `fh.E`.
The procedure yields fully correct results only if `fh.E` has been calculated at zero perturbation. The perturbation may additionally change
the order of levels, and there is no simple way of disentangling the order. The permutation is still useful in that case, but the results 
should not be taken too literally.
"""
function order_floquet_levels(fh::FloquetHamiltonian)
    E = similar(fh.E)
    n_levels = size(fh.E, 1) # number of Floquet levels; equivalently, number of levels of ℎ used for constructing ℋ
    for i in eachindex(fh.uh.φₓ)
        E_diag = [fh.uh.E[m, i] - fh.ν[m] * fh.ω/fh.s for m in 1:n_levels] # Floquet energies at zero perturbation
        invsort = sortperm(sortperm(E_diag, rev=true))  # inverse permutation, such that `sort(E_diag, rev=true)[invsort] == E_diag`
        E[:, i] .= fh.E[invsort, i]
    end
    return E
end

"""
Calculate Wannier vectors for the Floquet Hamiltonian `fh` using the quasienergy levels `targetlevels`.
"""
function compute_wanniers!(fh::FloquetHamiltonian; targetlevels::AbstractVector{<:Real})
    (;N, φₓ) = fh.uh

    n_w = length(targetlevels)
    E = Matrix{Float64}(undef, n_w, length(φₓ))
    pos = Matrix{Float64}(undef, n_w, length(φₓ))
    d = Array{ComplexF64, 3}(undef, n_w, n_w, length(φₓ))
    fh.uh.w = Wanniers(0, targetlevels, Int[], E, pos, d)
    X = Matrix{ComplexF64}(undef, n_w, n_w) # position operator
    
    n_levels = size(fh.E, 1)
    # matrices for storing intermediate results
    C = Matrix{ComplexF64}(undef, n_levels, n_levels)
    D = Matrix{ComplexF64}(undef, n_levels, n_levels)
    
    Ωt = π/5 # time moment at which to diagonalise the coordinate operator
    for i in eachindex(φₓ)
        # if pumping is time-only, then `C` must be calculated only at the first iteration, thereby using `c`'s at 𝜑ₓ = 0
        if fh.pumptype != :time || i == 1
            window = range(fh.minlevel, length=n_levels)
            mul!(C, @view(fh.uh.c[2:end, window, i])', @view(fh.uh.c[1:end-1, window, i]))
        end

        fh.pumptype != :space && (Ωt = π/5 - i/length(φₓ)*π/2)
        # `D` must be calculated at every phase: if pumping is temporal, `Ωt` depends on phase;
        # if pumping is spatial, `C` depends on phase (because `c`'s do)
        for m in 1:n_levels, m′ in 1:n_levels
            D[m′, m] = C[m′, m] * cis((fh.ν[m′+fh.minlevel-1] - fh.ν[m+fh.minlevel-1]) * Ωt)
        end

        X .= fh.b[:, targetlevels, i]' * D * fh.b[:, targetlevels, i]

        _, d[:, :, i], pos_complex = schur(X)
        pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
        sp = sortperm(pos_real)         # sort the eigenvalues
        pos[:, i] = pos_real[sp]
        @views Base.permutecols!!(d[:, :, i], sp) # sort the eigenvectors in the same way
        E[:, i] = transpose(fh.E[targetlevels, i]) * abs2.(d[:, :, i])
    end
end

"""
Construct Wannier functions at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `fh` are constructed.
In the process, energy eigenfunctions are also constructed.
Return `u, w`, where `w[ix, it, j, i]` = `j`th Wannier function at `whichphases[i]`th phase at `ix`th coordinate at `it`th time moment,
and `u` is an array of Floquet modes in the same format.
"""
function make_wannierfunctions(fh::FloquetHamiltonian, x::AbstractVector{<:Real}, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
    n_w = length(fh.uh.w.targetlevels)
    u = make_eigenfunctions(fh, x, Ωt, whichphases, fh.uh.w.targetlevels)
    w = Array{ComplexF64, 4}(undef, length(x), length(Ωt), n_w, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for j in 1:n_w
            w[:, :, j, i] = sum(fh.uh.w.d[k, j, iφ] * u[:, :, k, i] for k = 1:n_w)
        end
    end
    return u, w
end

"""
A type representing 2D (time+space) tight-binding Hamiltonian.
"""
mutable struct TBFloquetHamiltonian
    N::Int
    H::Array{ComplexF64, 3} # Hamiltonian matrix
    φₓ::Vector{Float64}
    pumptype::Symbol
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::Wanniers 
end

"""
Construct a `TBFloquetHamiltonian` object using the temporal band `targetband`.
`fh` must contain calculated periodic Wanniers; `pumptype` may or may not coincide with `fh.pumptype`.
"""
function TBFloquetHamiltonian(fh::FloquetHamiltonian; targetband::Integer, pumptype::Symbol)
    (;N, c, Vₗ, φₓ) = fh.uh
    (;s, λₗ, ν) = fh
    n_φₓ = length(φₓ)
    n_w = 4 * 2N # number of Wanniers
    H = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ) # TB Hamiltonian matrix
    
    iφ₀ = 1 # phase index at which to take the Wanniers -- any choice should work, but it is simpler to assume `iφ₀ = 1` below
    
    n_m = size(fh.E, 1) # number of levels of ℎ considered
    Ψ = Matrix{ComplexF64}(undef, n_m, n_m)
    
    d = @view fh.uh.w.d[:, :, iφ₀]
    bd = fh.b[:, range(n_w*(targetband-1) + 1, length=n_w), iφ₀] * d
    
    H₀ = d' * Diagonal(fh.E[range(n_w*(targetband-1) + 1, length=n_w), iφ₀]) * d

    for (iφ, φ) in enumerate(φₓ)
        for m in 1:n_m # `m` is the subband index of ℎ
            for m′ in 1:n_m
                Ψ[m′, m] = 0
                if pumptype != :time # if pumping is not time-only, account for the change of the spatial phase
                    if ν[m] == ν[m′]
                        @views Ψ[m′, m] += Vₗ / 4 * ( (c[1+N:end, m′, iφ₀] ⋅ c[1:end-N, m, iφ₀]) * (cis(+2φ) - 1) +
                                                      (c[1:end-N, m′, iφ₀] ⋅ c[1+N:end, m, iφ₀]) * (cis(-2φ) - 1) )
                    end
                end
                if pumptype != :space # if pumping is not space-only, account for the change of the temporal phase
                    if ν[m] == ν[m′] + s
                        e = cis(+2φ) # we assume 𝜑ₜ = 2𝜑ₓ, hence the two
                    elseif ν[m] == ν[m′] - s
                        e = cis(-2φ)
                    else
                        continue
                    end
                    @views Ψ[m′, m] += λₗ/8 * (e - 1) * ( (c[1+2N:end, m′, iφ₀] ⋅ c[1:end-2N, m, iφ₀]) +
                                                          (c[1:end-2N, m′, iφ₀] ⋅ c[1+2N:end, m, iφ₀]) )
                end
            end
        end
        H[:, :, iφ] = H₀ + bd' * Ψ * bd
    end

    E = Matrix{Float64}(undef, n_w, n_φₓ)
    cc = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ)
    TBFloquetHamiltonian(N, H, φₓ, pumptype, E, cc, Wanniers())
end

"Diagonalise the TB Floquet Hamiltonian `tbh` at each phase."
function diagonalise!(tbh::TBFloquetHamiltonian)
    for iφ in eachindex(tbh.φₓ)
        tbh.E[:, iφ], tbh.c[:, :, iφ] = eigen(Hermitian(tbh.H[:, :, iφ]))
    end
end

end