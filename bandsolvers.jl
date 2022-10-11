module Bandsolvers

import BandedMatrices as BM
using SparseArrays: sparse
using KrylovKit: eigsolve
using LinearAlgebra: eigen, eigvals, schur, ⋅, diagm, diagind

"A type representing the spatial Wannier functions."
mutable struct SpatialWanniers
    minlevel::Int # number of the first energy level to use for constructing wanniers
    E_lo::Vector{Vector{Float64}} # `E[j][i]` = mean energy of `i`th wannier at `j`th phase; use nested arrays because in non-periodic case
    E_hi::Vector{Vector{Float64}} #               the number of wanniers in each subband varies depending on the edge state branch position
    pos_lo::Vector{Vector{Float64}} # `pos[j][i]` = position eigenvalue of `i`th wannier at `j`th phase
    pos_hi::Vector{Vector{Float64}}
    d_lo::Vector{Matrix{ComplexF64}} # `d[j][:, i]` = 𝑖th wannier vector at `j`th phase
    d_hi::Vector{Matrix{ComplexF64}}
end

"""
A type representing the unperturbed Hamiltonian
    ℎ = 𝑝/2𝑀 + 𝑔ₗcos²(2𝑥) + 𝑉ₗcos²(𝑥 + 𝜑ₓ).
"""
mutable struct UnperturbedHamiltonian
    N::Int # number of lattice cells
    M::Float64
    l::Int
    gₗ::Float64
    Vₗ::Float64
    isperiodic::Bool
    phases::Vector{Float64}   # values of the spatial adiabatic phases 𝜑ₓ
    maxlevel::Int # highest level number to consider
    bandsizes::Tuple{Int, Int} # = (number of levels in the first band, number of levels in the second band)
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, `i = 1:maxlevel`
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase, `i = 1:maxlevel; j = 1:2maxlevel+1`
    w::SpatialWanniers
end

"""
Construct an `UnperturbedHamiltonian` object. `maxband` is the highest energy band number to consider.
Each band is assumed to consist of 2 subbands containing `2n_cells` levels in total.
"""
function UnperturbedHamiltonian(n_cells::Integer; M::Real, gₗ::Real, Vₗ::Real, maxband::Integer, isperiodic::Bool, phases::AbstractVector{<:Real},
                                l::Union{Nothing, Integer}=nothing)
    bandsizes = (2n_cells - 1, 2n_cells + 1)

    # convert max band number to level number
    if isperiodic
        maxlevel = maxband * 2n_cells
    else
        maxlevel = (maxband-1) ÷ 2 * 4n_cells + (isodd(maxband) ? bandsizes[1] : sum(bandsizes))
    end

    E = Matrix{Float64}(undef, maxlevel, length(phases))
    c = Array{ComplexF64,3}(undef, 2maxlevel+1, maxlevel, length(phases))

    E_lo = [Float64[] for _ in eachindex(phases)]
    E_hi = [Float64[] for _ in eachindex(phases)]
    pos_lo = [Float64[] for _ in eachindex(phases)]
    pos_hi = [Float64[] for _ in eachindex(phases)]
    d_lo = [ComplexF64[;;] for _ in eachindex(phases)]
    d_hi = [ComplexF64[;;] for _ in eachindex(phases)]
    w = SpatialWanniers(0, E_lo, E_hi, pos_lo, pos_hi, d_lo, d_hi)

    UnperturbedHamiltonian(Int(n_cells), Float64(M), (l === nothing ? 1 : l), Float64(gₗ), Float64(Vₗ), isperiodic,
                           collect(Float64, phases), maxlevel, bandsizes, E, c, w)
end

"Diagonalise the unperturbed Hamiltonian at each phase."
function diagonalise!(uh::UnperturbedHamiltonian)
    (;N, M, gₗ, Vₗ, maxlevel) = uh
    sortby = M > 0 ? (+) : (-) # eigenvalue sorting; for 𝑀 < 0 we use descending sorting
    if uh.isperiodic
        h = diagm(0 => ComplexF64[(2j/N)^2 / 2M + (gₗ + Vₗ)/2 for j = -maxlevel:maxlevel])
        h[diagind(h, -2N)] .= h[diagind(h, 2N)] .= gₗ/4
        for (i, ϕ) in enumerate(uh.phases)
            h[diagind(h, -N)] .= Vₗ/4 * cis(+2ϕ)
            h[diagind(h, +N)] .= Vₗ/4 * cis(-2ϕ)
            f = eigen(h; sortby)
            uh.E[:, i] = f.values[1:maxlevel]
            uh.c[:, :, i] = f.vectors[:, 1:maxlevel]
        end
    else
        X(j′, j) = 16N*j*j′ / (π*((j-j′)^2-(2N)^2)*((j+j′)^2-(2N)^2))
        n_j = 2maxlevel + 1
        h = zeros(n_j, n_j)
        for (i, ϕ) in enumerate(uh.phases)
            for j in 1:n_j
                for j′ in j:n_j
                    val = 0.0
                    if isodd(j′ + j)
                        val += Vₗ/2 * X(j′, j) * sin(2ϕ)
                    else
                        # check diagonals "\"
                        if j′ == j
                            val += (gₗ + Vₗ)/2 + (j / N)^2 / 2M
                        elseif j′ == j - 2N || j′ == j + 2N
                            val += Vₗ * cos(2ϕ) / 4
                        elseif j′ == j - 4N || j′ == j + 4N
                            val += gₗ/4
                        end
                        # check anti-diagonals "/"
                        if j′ == -j - 2N || j′ == -j + 2N
                            val += -Vₗ * cos(2ϕ) / 4
                        elseif j′ == -j - 4N || j′ == -j + 4N
                            val += -gₗ/4
                        end
                    end
                    h[j′, j] = h[j, j′] = val # push the element to the conjugate positions
                end
            end
            f = eigen(h; sortby)
            uh.E[:, i] = f.values[1:maxlevel]
            uh.c[:, :, i] = f.vectors[:, 1:maxlevel]
        end
    end
end

"""
Calculate Wannier vectors for the unperturbed Hamiltonian using the energy eigenstates in the band number `targetband`.
Note that if mass is negative (`uh.M < 0`), then `uh.w.E_lo`, `uh.w.pos_lo`, and `uh.w.pos_lo` will refer to the band whose energy is higher.
"""
function compute_wanniers!(uh::UnperturbedHamiltonian; targetband::Integer)
    N = uh.N

    minlevel = (targetband-1) * 2N + 1
    uh.w.minlevel = minlevel # save this because it's needed in `make_wavefunction` when constructing coordinate space Wannier functions

    if uh.isperiodic
        X = Matrix{ComplexF64}(undef, N, N) # position operator
        
        for i in eachindex(uh.phases)
            # Lower band
            for n in 1:N
                for n′ in 1:N
                    X[n′, n] = sum(uh.c[j+1, minlevel+n′-1, i]' * uh.c[j, minlevel+n-1, i] for j = 1:size(uh.c, 1)-1)
                end
            end
            # `eigen` does not guarantee orthogonality of eigenvectors in case of degeneracies for `X` unitary, so use `schur`
            # (although a degeneracy of coordinates eigenvalues is unlikely here)
            _, uh.w.d_lo[i], pos_complex = schur(X)
            pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
            sp = sortperm(pos_real)                 # sort the eigenvalues
            uh.w.pos_lo[i] = pos_real[sp]
            Base.permutecols!!(uh.w.d_lo[i], sp)    # sort the eigenvectors in the same way
            uh.w.E_lo[i] = [abs2.(dˣ) ⋅ uh.E[minlevel:minlevel+N-1, i] for dˣ in eachcol(uh.w.d_lo[i])]

            # Higher band
            for n in 1:N
                for n′ in 1:N
                    X[n′, n] = sum(uh.c[j+1, minlevel+N+n′-1, i]' * uh.c[j, minlevel+N+n-1, i] for j = 1:size(uh.c, 1)-1)
                end
            end
            _, uh.w.d_hi[i], pos_complex = schur(X)
            pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
            sp = sortperm(pos_real)                 # sort the eigenvalues
            uh.w.pos_hi[i] = pos_real[sp]
            Base.permutecols!!(uh.w.d_hi[i], sp)    # sort the eigenvectors in the same way
            uh.w.E_hi[i] = [abs2.(dˣ) ⋅ uh.E[minlevel+N:minlevel+2N-1, i] for dˣ in eachcol(uh.w.d_hi[i])]
        end
    else
        n_w = isodd(targetband) ? 2N-1 : 2N+1 # total number of wanniers to construct; this is the number of levels in the target band
        X_less = zeros(n_w÷2, n_w÷2) # position operator for a subband which does not contain the edge state branch
        X_more = zeros(n_w÷2+1, n_w÷2+1) # position operator for a subband which contains the edge state branch
        
        n_j = size(uh.c, 1)

        for i in eachindex(uh.phases)
            up = uh.E[minlevel+n_w÷2, i] > (uh.E[minlevel+n_w÷2-1, i] + uh.E[minlevel+n_w÷2+1, i])/2 # true if the edge state branch is above the mean value
            
            # Lower band
            X = up ? X_less : X_more # bind the position operator to a matrix of the appropriate size
            n_lo = n_w÷2 + !up # number of levels in the lower subband
            for n in 1:n_lo
                for n′ in n:n_lo
                    X[n′, n] = X[n, n′] = N*sum(uh.c[j, minlevel+n-1, i] * (π/2 * uh.c[j, minlevel+n′-1, i] - 8/π * sum(uh.c[j′, minlevel+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                            for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
                end
            end
            uh.w.pos_lo[i], uh.w.d_lo[i] = eigen(X)
            uh.w.E_lo[i] = [dˣ.^2 ⋅ uh.E[minlevel:minlevel+n_lo-1, i] for dˣ in eachcol(uh.w.d_lo[i])]

            # Higher band
            X = up ? X_more : X_less
            n_hi = n_w - n_lo
            for n in 1:n_hi
                for n′ in n:n_hi
                    X[n′, n] = X[n, n′] = N*sum(uh.c[j, minlevel+n_lo+n-1, i] * (π/2 * uh.c[j, minlevel+n_lo+n′-1, i] - 8/π * sum(uh.c[j′, minlevel+n_lo+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                                 for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
                end
            end
            uh.w.pos_hi[i], uh.w.d_hi[i] = eigen(X)
            uh.w.E_hi[i] = [dˣ.^2 ⋅ uh.E[minlevel+n_lo:minlevel+n_w-1, i] for dˣ in eachcol(uh.w.d_hi[i])]
        end
    end
end

"""
Construct energy eigenfunctions `ψ` at coordinates in `x` for each eigenstate number in `whichstates` at each phase number in `whichphases`.
Return `ψ`, where `ψ[:, j, i]` = `j`th wavefunction at `i`th phase.
"""
function make_eigenfunctions(uh::UnperturbedHamiltonian, x::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer}, whichstates::AbstractVector{<:Integer})
    ψ = Array{ComplexF64,3}(undef, length(x), length(whichstates), length(whichphases))
    make_state = uh.isperiodic ? make_exp_state : make_sin_state
    for (i, iϕ) in enumerate(whichphases)
        for (j, js) in enumerate(whichstates)
            ψ[:, j, i] = make_state(x, uh.c[:, js, iϕ]; N=uh.N)
        end
    end
    return ψ
end

"""
Construct Wannier functions `w_lo` and `w_hi` at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `uh` are constructed.
Return `(w_lo, w_hi)`, where `w_xx[i][j]` = `j`th Wannier function at `i`th phase.
"""
function make_wannierfunctions(uh::UnperturbedHamiltonian, x::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
    w_lo = [Vector{Vector{ComplexF64}}() for _ in eachindex(whichphases)]
    w_hi = [Vector{Vector{ComplexF64}}() for _ in eachindex(whichphases)]
    make_state = uh.isperiodic ? make_exp_state : make_sin_state
    for (i, iϕ) in enumerate(whichphases)
        n_lo = size(uh.w.d_lo[iϕ], 1)
        w_lo[i] = [Vector{ComplexF64}(undef, length(x)) for _ in 1:n_lo]
        for j in 1:n_lo
            w_lo[i][j] = sum(uh.w.d_lo[iϕ][k, j] * make_state(x, uh.c[:, uh.w.minlevel+k-1, iϕ]; N=uh.N) for k = 1:n_lo)
        end
        n_hi = size(uh.w.d_hi[iϕ], 1)
        w_hi[i] = [Vector{ComplexF64}(undef, length(x)) for _ in 1:n_hi]
        for j in 1:n_hi
            w_hi[i][j] = sum(uh.w.d_hi[iϕ][k, j] * make_state(x, uh.c[:, uh.w.minlevel+n_lo+k-1, iϕ]; N=uh.N) for k = 1:n_hi)
        end
    end
    return w_lo, w_hi
end

"""
Reshape the Wannier function `w` returned by [`make_wannierfunctions`](@ref) into a matrix for a fixed Wannier state number `n`.
Return `w_map`, where `w_map[:, i]` = `n`th Wannier function at `i`th phase.
"""
function make_wanniermap(w::Vector{Vector{Vector{T}}}, n::Integer) where T <: Number
    w_map = Matrix{T}(undef, length(w[1][1]), length(w))
    for i in eachindex(w)
        w_map[:, i] = w[i][n]
    end
    return w_map
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
A type representing the Floquet Hamiltonian
    ℋ = ℎ - i∂ₜ + λₛsin²(2𝑥)cos(2𝜔𝑡) + λₗcos²(2𝑥)cos(𝜔𝑡 + 𝜑ₜ),
where ℎ is the unperturbed Hamiltonian represented by [`UnperturbedHamiltonian`](@ref).
Type of pumping is controlled via `pumptype`: `:time` for temporal, `:space` for spatial, or anything else for simultaneous space-time pumping.
"""
mutable struct FloquetHamiltonian
    uh::UnperturbedHamiltonian
    s::Int
    λₛ::Float64
    λₗ::Float64
    ω::Float64
    pumptype::Symbol
    minlevel::Int # lowest level number of ℎ to use when constructing Floquet Hamiltonian
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue (Floquet quasienergy) at `j`th phase, `i = 1:maxlevel`
    b::Array{ComplexF64, 3} # `b[:, i, j]` = `i`th eigenvector at `j`th phase, `i = 1:maxlevel; j = 1:2maxlevel+1`
end

function FloquetHamiltonian(uh::UnperturbedHamiltonian; s::Integer, λₛ::Real, λₗ::Real, ω::Real, pumptype::Symbol, minband::Integer)
    N = uh.N

    # convert min band number to level number
    if uh.isperiodic
        minlevel = (minband - 1) * 2N + 1
    else
        minlevel = (minband - 1) ÷ 2 * 4N + (isodd(minband) ? 1 : uh.bandsizes[1] + 1)
    end

    # if iseven(minlevel) # swap `bs1` and `bs2` so that they correspond to actual band sizes
    #     bs1, bs2 = bs2, bs1
    # end

    n_levels = uh.maxlevel - minlevel + 1
    # ν = Vector{Int}(undef, Δn)
    # # FIll `ν`: [1 (`bs1` times), 2 (`bs2` times), 3 (`bs1` times), 4 (`bs2` times), ...]
    # number = 1
    # g = bs1 + bs2
    # for i in 0:Δn÷g-1
    #     ν[g*i+1:g*i+bs1] .= number
    #     number += 1
    #     ν[g*i+bs1+1:g*i+g] .= number
    #     number += 1
    # end
    # ν[Δn - Δn%g + 1:end] .= number

    # pattern = [fill(bs1, bs1); fill(bs2, bs2)]
    # G = repeat(pattern, Δn÷g) # a pattern which e.g. for `n == 2` looks like [3, 3, 3, 5, 5, 5, 5, 3, 3, 3, 5, 5, 5, 5, ...]
    # Δn % g != 0 && append!(G, fill(bs1, bs1))
    
    E = Matrix{Float64}(undef, n_levels, length(uh.phases))
    b = Array{ComplexF64,3}(undef, n_levels, n_levels, length(uh.phases))
    
    FloquetHamiltonian(uh, Int(s), Float64(λₛ), Float64(λₗ), Float64(ω), pumptype, Int(minlevel), E, b)
end

"""
Diagonalise Floquet Hamiltonian for a periodic system and calculate Wannier centres.
    `n_min` - lowest band number of spatial Hamiltonian to use when constructing Floquet Hamiltonian
    `n_max` - highest band number of spatial Hamiltonian to consider
    `n_target` - number of Floquet band (counting from the highest) to use for Wannier calculations
    `coords` - x's for wavefunctions
    `ωts` - time moments for wavefunctions
    `mix_time_cells` - if set to `false`, only the spatial levels of the highest temporal level will be used when constructing Wanniers
"""
function diagonalise!(fh::FloquetHamiltonian)
    # n_target_min = (n_target-1) * 4N + 1

    # n_j = n_max * 2N # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian

    # h = zeros(ComplexF64, 2n_j + 1, 2n_j + 1)
    # h = diagm(0 => ComplexF64[(2j/N)^2 + (gₗ + Vₗ)/2 for j = -n_j:n_j])
    # h[diagind(h, -2N)] .= h[diagind(h, 2N)] .= gₗ/4

    # n_w = mix_time_cells ? s*N : N # number of Wannier functions to construct
    # pos_lo = Matrix{Float64}(undef, n_w, length(phases)) # position eigenvalues (Wannier centres) of the lower spatial levels
    # pos_hi = Matrix{Float64}(undef, n_w, length(phases)) # position eigenvalues (Wannier centres) of the higher spatial levels
    # ε_lo = Matrix{Float64}(undef, n_w, length(phases)) # energies of Wannier states of the lower spatial levels
    # ε_hi = Matrix{Float64}(undef, n_w, length(phases)) # energies of Wannier states of the higher spatial levels
    # wf_lo = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the lower spatial levels
    # wf_hi = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the higher spatial levels
    # window_lo = Int[] # window of Floquet eigenstates which will be used to construct Wannier functions out of the lower spatial levels
    # window_hi = Int[] # window of Floquet eigenstates which will be used to construct Wannier functions out of the higher spatial levels
    # for i in 0:n_w÷N - 1
    #     append!(window_hi, n_target_min+2i*N:n_target_min+2i*N + N - 1)
    #     append!(window_lo, n_target_min+(2i+1)*N:n_target_min+(2i+1)*N + N - 1)
    # end

    # u_lo = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the lower spatial levels
    # u_hi = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the higher spatial levels

    # make views
    (;N, phases, E, c) = fh.uh
    n_j = size(c, 1)
    (;s, ω, λₛ, λₗ, pumptype) = fh

    n_levels = fh.uh.maxlevel - fh.minlevel + 1 # number of levels of spatial Hamiltonian to use for constructing Floquet Hamiltonian
    ν(m) = ceil(Int, m/2N)

    # ψ = Matrix{ComplexF64}(undef, length(coords), n_levels) # eigenvectors of ℎ in 𝑥-representation
    # cc = Matrix{ComplexF64}(undef, n_levels, n_levels) # matrix of products of `c`'s that will be needed multiple times
    # ccc = Matrix{ComplexF64}(undef, n_levels, n_levels) # products `cc`'s and cis
    
    # x = Matrix{ComplexF64}(undef, n_w, n_w) # position operator

    H = zeros(ComplexF64, n_levels, n_levels) # Floquet Hamiltonian matrix

    for (i, ϕ) in enumerate(fh.uh.phases)
            # # construct coordinate representation of eigenfunctions and compute products of `c`'s that will be needed multiple times
            # for m in 1:n_levels
            #     ψ[:, m] = make_exp_state(coords, c[:, m]; n=N)
            #     for m′ in 1:n_levels
            #         cc[m′, m] = sum(c[j+1, m′]' * c[j, m] for j = 1:2n_j)
            #     end
            # end
            # if pumptype == :time
            #     for p in 2:length(phases) # copy the calculated first column of `ϵ` to all other columns for consistency
            #         ϵ[:, p] = ϵ[:, 1]
            #     end
            # end

        # `m` and `m′` number the levels of the unperturbed Hamiltonian
        # `e` and `e′` number the elements of the FloquetHamiltonian
        for m in fh.minlevel:fh.uh.maxlevel
            e = m - fh.minlevel + 1

            # for time-only pumping, always take the eigenenergies at the first phase, corresponding to 𝜑ₓ = 0
            p = (pumptype == :time ? 1 : i)
            H[e, e] = E[m, p] - ν(m)*ω/s

            # place the elements of the long lattice
            for g in 1:2N
                m′ = 2N*(s + ν(m) - 1) + g
                e′ = m′ - fh.minlevel + 1
                e′ > n_levels && break
                if pumptype != :time || i == 1 # if pumping is time-only, this may be calculated only once, at `i` = 1
                    j_sum = sum( (                 2c[j, m′, i] + c[j+2N, m′, i])' * c[j, m, i] for j = 1:2N ) +
                            sum( (c[j-2N, m′, i] + 2c[j, m′, i] + c[j+2N, m′, i])' * c[j, m, i] for j = 2N+1:n_j-2N ) + 
                            sum( (c[j-2N, m′, i] + 2c[j, m′, i]                 )' * c[j, m, i] for j = n_j-2N+1:n_j )
                    # if pumping is space-time, then also multiply by cis(-𝜑ₜ). `ϕ` runs over the spatial phases 𝜑ₓ,
                    H[e′, e] = (pumptype == :space ? λₗ/8 * j_sum : λₗ/8 * j_sum * cis(-2ϕ)) # and we assume the pumping protocol 𝜑ₜ = 2𝜑ₓ
                elseif pumptype == :time 
                    H[e′, e] *= cis(-2(phases[2]-phases[1]))
                end
                H[e, e′] = H[e′, e]'
            end
            
            # place the elements of the short lattice
            for g in 1:2N
                m′ = 2N*(2s + ν(m) - 1) + g
                e′ = m′ - fh.minlevel + 1
                e′ > n_levels && break
                if pumptype != :time || i == 1 # if pumping is time-only, this may be calculated only once, at `i` = 1
                    j_sum = sum( (                  2c[j, m′, i] - c[j+2N, m′, i])' * c[j, m, i] for j = 1:2N ) +
                            sum( (-c[j-2N, m′, i] + 2c[j, m′, i] - c[j+2N, m′, i])' * c[j, m, i] for j = 2N+1:n_j-2N ) + 
                            sum( (-c[j-2N, m′, i] + 2c[j, m′, i]                 )' * c[j, m, i] for j = n_j-2N+1:n_j )
                    H[e′, e] = λₛ/8 * j_sum
                end
                H[e, e′] = H[e′, e]'
            end
        end
        fh.E[:, i], fh.b[:, :, i] = eigen(H, sortby=-)

        # ### Wannier centres
        
        # t = (pumptype == :space ? π/5 : π/5 - z/length(phases)*π/2) # time moment at which to diagonalise the coordinate operator
        # for m in 1:n_levels, m′ in 1:n_levels
        #     ccc[m′, m] = cc[m′, m] * cis((ν(m′)-ν(m))*t)
        # end

        # # Higher band
        # # the loop below runs faster if we make a copy rather than a view of `f.vectors`; 
        # # both approaches are ~6 times faster compared to iterating directly over `f.vectors`
        # b .= f.vectors[:, window_hi]
        # for n in 1:n_w, n′ in 1:n_w
        #     x[n′, n] = sum(b[m, n] * sum(b[m′, n′]' * ccc[m′, m] for m′ in 1:n_levels) for m in 1:n_levels)
        # end
        # _, d, pos_complex = schur(x)
        # pos_real = (angle.(pos_complex) .+ π) / 2π * N*π # take angle, convert from (-π, π) to (0, 2π), and map to the interval (0, Nπ)
        # sp = sortperm(pos_real)
        # pos_hi[:, z] = pos_real[sp]   # sort positions
        # Base.permutecols!!(d, sp)     # sort eigenvectors in the same way
        # ε_hi[:, z] = [abs2.(dˣ) ⋅ E[window_hi, z] for dˣ in eachcol(d)]
        # for (t, ωt) in enumerate(ωts)
        #     for X in 1:n_w
        #         wf_hi[:, X, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * sum(d[l, X] * b[m, l] for l = 1:n_w) for m in 1:n_levels))
        #     end
        # end
        # # for (t, ωt) in enumerate(ωts)
        # #     for l in 1:n_w
        # #         u_hi[:, l, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * b[m, l] for m in 1:n_levels))
        # #     end
        # # end

        # # Lower band
        # b .= f.vectors[:, window_lo]
        # for n in 1:n_w, n′ in 1:n_w
        #     x[n′, n] = sum(b[m, n] * sum(b[m′, n′]' * ccc[m′, m] for m′ in 1:n_levels) for m in 1:n_levels)
        # end
        # _, d, pos_complex = schur(x)
        # pos_real = (angle.(pos_complex) .+ π) / 2π * N*π # take angle, convert from (-π, π) to (0, 2π), and map to the interval (0, Nπ)
        # sp = sortperm(pos_real)
        # pos_lo[:, z] = pos_real[sp]   # sort positions
        # Base.permutecols!!(d, sp)     # sort eigenvectors in the same way
        # ε_lo[:, z] = [abs2.(dˣ) ⋅ E[window_lo, z] for dˣ in eachcol(d)]
        # for (t, ωt) in enumerate(ωts)
        #     for X in 1:n_w
        #         wf_lo[:, X, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * sum(d[l, X] * b[m, l] for l = 1:n_w) for m in 1:n_levels))
        #     end
        # end
    end
end

"""
Construct Floquet modes at coordinates in `x` and time moments in `Ωt` for each tate number in `whichstates` at each phase number in `whichphases`.
Return `u`, where `u[ix, it, j, i]` = `j`th wavefunction at `i`th phase at `ix`th coordinate at `it`th time moment.
"""
function make_eigenfunctions(fh::FloquetHamiltonian, x::AbstractVector{<:Real}, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer},
                             whichstates::AbstractVector{<:Integer})
    u = Array{ComplexF64,4}(undef, length(x), length(Ωt), length(whichstates), length(whichphases))
    n_levels = size(fh.E, 1) # number of Floquet levels; equivalently, number of levels of ℎ used for constructing ℋ
    ψ = Array{ComplexF64,2}(undef, length(x), n_levels) # for storing eigenfunctions of ℎ, which are mixed during construction of `u`
    ν(m) = ceil(Int, m/2fh.uh.N)
    make_state = fh.uh.isperiodic ? make_exp_state : make_sin_state
    for (i, iϕ) in enumerate(whichphases)
        for m in 1:n_levels
            ψ[:, m] = make_state(x, fh.uh.c[:, fh.minlevel+m-1, iϕ]; N=fh.uh.N)
        end
        for (j, js) in enumerate(whichstates)
            for (it, t) in enumerate(Ωt)
                u[:, it, j, i] = sum(cis(-ν(m)*t) * ψ[:, m] * fh.b[m, js, iϕ] for m in 1:n_levels)
            end
        end
    end
    return u
end

end

"""
Calculate energy bands of the Floquet Hamiltonian (S20) sweeping over the adiabatic `phases` φₓ. It is assumed that 2φₜ = φₓ.
Energy levels of the unperturbed Hamiltonian ℎₖ from `2n_min` to `2n_max` will be used for constructing the Floquet Hamiltonian.
The values `n_min` to `n_max` thus correspond to the energy level numbers of a single well.
Return a tuple of a matrix `ϵₖ` of `4Δn` bands of ℎₖ and a matrix `Eₖ` of `Δn` bands of 𝐻ₖ, where `Δn = n_max - n_min + 1`.
In the returned matrices, columns enumerate the adiabatic phases, while rows enumerate eigenvalues.
In `Eₖ`, rows `1:Δn` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
In `Eₖ`, rows `Δn:end` store the eigenvalues corresponding to the boundary of BZ, in our case Vₗcos²(x+φₓ) leads to 𝑘 = 2/2 = 1.
The structure of `ϵₖ` is the same, but with `2Δn` instead of `Δn`.
Type of pumping is controlled via `pumptype`: `:time` for temporal, `:space` for spatial, or anything else for simultaneous space-time pumping.
"""
function compute_floquet_bands(; n_min::Integer, n_max::Integer, phases::AbstractVector{<:Real}, s::Integer, l::Real, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
    n_j = 2n_max # number of indices 𝑗 to use for constructing ℎₖ (its size will be (2n_j+1)×(2n_j+1)). `2n_max` is a safe value, but it could be less.
    Δn = n_max - n_min + 1

    hₖ = BM.BandedMatrix(BM.Zeros{ComplexF64}(2n_j+1, 2n_j+1), (2l, 2l))
    # fill the off-diagonals with binomial numbers; the diagonal is treated in the `k` loop
    for n in 1:l
        hₖ[BM.band(2n)] .= hₖ[BM.band(-2n)] .= gₗ / 4^l * binomial(2l, l-n)
    end
    
    # Eigenvalues of ℎₖ (eigenenergies of the unperturbed Hamiltonian).
    # We should store `2Δn` of them because each of the `Δn` levels are almost degenerate. To account for the two values of 𝑘, we use `4Δn`.
    ϵₖ = Matrix{Float64}(undef, 4Δn, length(phases))
    cₖ = [Vector{ComplexF64}(undef, 2n_j+1) for _ in 1:2Δn]  # eigenvectors of ℎₖ, we will save `2Δn` of them (only for a single 𝑘), and each will have `2n_j+1` components
    
    Eₖ = Matrix{Float64}(undef, 4Δn, length(phases)) # eigenvalues of 𝐻ₖ (Floquet quasi-energies) that will be saved; size is twice `Δn` for the two values of 𝑘
    Hₖ_dim = 2Δn # dimension of the constructed 𝐻ₖ matrix (twice larger than the number of requested quasi-energies)
    n_Hₖ_nonzeros = 9Hₖ_dim - 24s # number of non-zero elements in 𝐻ₖ
    Hₖ_rows = Vector{Int}(undef, n_Hₖ_nonzeros)
    Hₖ_cols = Vector{Int}(undef, n_Hₖ_nonzeros)
    Hₖ_vals = Vector{ComplexF64}(undef, n_Hₖ_nonzeros)
    for k in [0, 1] # iterate over the centre of BZ and then the boundary
        hₖ[BM.band(0)] .= [(2j + k)^2 + Vₗ/2 + gₗ / 4^l * binomial(2l, l) for j = -n_j:n_j]
        # `a_*` and `b_*` control where to place the eigenvalues of 𝐻ₖ and ℎₖ depedning on `k`; see function docstring
        a_hₖ = (k > 0)*2Δn + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_hₖ = a_hₖ+2Δn - 1
        a_Hₖ = (k > 0)*2Δn + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_Hₖ = a_Hₖ+2Δn - 1
        for (z, ϕ) in enumerate(phases)
            if pumptype != :time || z == 1 # If pupming is not time-only, ℎₖ has to be diagonalised on each iteration. If it's time-only, then we diagonalise only once, at `z == 1`.
                hₖ[BM.band(-1)] .= Vₗ/4 * cis(2ϕ)
                hₖ[BM.band(1)]  .= Vₗ/4 * cis(-2ϕ)
                vals, vecs, info = eigsolve(hₖ, 2n_max, :SR; krylovdim=2n_j+1)
                if info.converged < 2n_max
                    @warn "Only $(info.converged) eigenvalues out of $(2n_max) converged when diagonalising ℎₖ. "*
                          "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
                end
                # save only energies and states for levels from `2n_min-1` to `2n_max`
                ϵₖ[a_hₖ:b_hₖ, z] = vals[2n_min-1:2n_max]
                cₖ .= vecs[2n_min-1:2n_max]
                if pumptype == :time
                    for p in 2:length(phases) # copy the calculated first column of `ϵₖ` to all other columns for consistency
                        ϵₖ[a_hₖ:b_hₖ, p] = ϵₖ[a_hₖ:b_hₖ, 1]
                    end
                end
            end

            # Construct 𝐻ₖ
            p = 1 # a counter for placing elements to the vectors `Hₖ_*`
            for m in 1:Hₖ_dim
                # place the diagonal element (S25)
                Hₖ_rows[p] = Hₖ_cols[p] = m
                Hₖ_vals[p] = ϵₖ[m+a_hₖ-1, z] - ceil((m+2n_min-2)/2)*ω/s
                p += 1

                # place the elements of the long lattice (S26)
                for i in 1:2
                    m′ = 2s + 2(ceil(Int, m/2)-1) + i
                    m′ > Hₖ_dim && break
                    Hₖ_rows[p] = m′
                    Hₖ_cols[p] = m
                    if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                        j_sum = sum( (                cₖ[m′][j]/2 + cₖ[m′][j+2]/4)' * cₖ[m][j] for j = 1:2 ) +
                                sum( (cₖ[m′][j-2]/4 + cₖ[m′][j]/2 + cₖ[m′][j+2]/4)' * cₖ[m][j] for j = 3:2n_j-1 ) + 
                                sum( (cₖ[m′][j-2]/4 + cₖ[m′][j]/2                )' * cₖ[m][j] for j = 2n_j:2n_j+1)
                        Hₖ_vals[p] = (pumptype == :space ? λₗ/2 * j_sum : λₗ/2 * j_sum * cis(-2ϕ)) # a check for space or space-time pumping
                    elseif pumptype == :time 
                        Hₖ_vals[p] *= cis(-2(phases[2]-phases[1]))
                    end
                    p += 1
                    # place the conjugate element
                    Hₖ_rows[p] = m
                    Hₖ_cols[p] = m′
                    Hₖ_vals[p] = Hₖ_vals[p-1]'
                    p += 1
                end
                
                # place the elements of the short lattice (S29)
                for i in 1:2
                    m′ = 4s + 2(ceil(Int, m/2)-1) + i
                    m′ > Hₖ_dim && break
                    Hₖ_rows[p] = m′
                    Hₖ_cols[p] = m
                    if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                        j_sum = sum( (                 cₖ[m′][j]/2 - cₖ[m′][j+2]/4)' * cₖ[m][j] for j = 1:2 ) +
                                sum( (-cₖ[m′][j-2]/4 + cₖ[m′][j]/2 - cₖ[m′][j+2]/4)' * cₖ[m][j] for j = 3:2n_j-1 ) + 
                                sum( (-cₖ[m′][j-2]/4 + cₖ[m′][j]/2                )' * cₖ[m][j] for j = 2n_j:2n_j+1)
                        Hₖ_vals[p] = λₛ/2 * j_sum
                    end
                    p += 1
                    # place the conjugate element
                    Hₖ_rows[p] = m
                    Hₖ_cols[p] = m′
                    Hₖ_vals[p] = Hₖ_vals[p-1]'
                    p += 1
                end
            end
            Hₖ = sparse(Hₖ_rows, Hₖ_cols, Hₖ_vals)
            vals, _, info = eigsolve(Hₖ, 2Δn, :LR; krylovdim=Hₖ_dim)
            if info.converged < 2Δn
                @warn "Only $(info.converged) eigenvalues out of $(2Δn) converged when diagonalising 𝐻ₖ. "*
                      "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
            end
            Eₖ[a_Hₖ:b_Hₖ, z] .= vals[1:2Δn]
        end
    end
    return ϵₖ, Eₖ
end

"""
Permute Floquet energy levels (i.e. rows) contained in `E` so that they are stored in the same order as the eigenenergies `e` of
the spatial Hamiltonian (i.e. in the order of increasing 𝑚). Repeat this for every phase (i.e. column of `E`) and for both halves
of rows of `E` (first half is 𝑘 = 0, and second is 𝑘 = 1).
It is assumed that the bands in `E` are initially stored in energy-descending order, as obtained during diagonalisation. To perfrorm
the sorting, we first calculate `e - ν(m)` which is the diagonal of the Floquet Hamiltonian. If there is no perturbation, then these
are the Floquet quasienergies. Then, we sort them in descending order (as if we diagonalised the Hamiltonian) and find the permutation
that would undo this sorting. This permutation is applied to `E`.
The procedure yields fully correct results only if `E` has been calculated at zero perturbation. The perturbation may additionally change
the order of levels, and there is no simple way to disentangle the order. The permutation is still useful in that case, but the results 
should not be taken too literally.
"""
function permute_floquet_bands!(E::AbstractMatrix{<:Float64}, e::AbstractMatrix{<:Float64}, n_min::Integer, ω::Real, s::Integer)
    ν(m) = ceil((m+2n_min-2)/2) * ω/s
    n_energies = size(e, 1) ÷ 2 
    n_phases = size(e, 2)
    for p in 1:n_phases
        for k in [0, 1] # iterate over the centre of BZ and then the boundary
            offset = (k > 0)*n_energies # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
            e_diag = [e[m+offset, p] - ν(m) for m in 1:n_energies] # Floquet energies at zero perturbation
            invsort = sortperm(sortperm(e_diag, rev=true)) .+ offset # inverse permutation, such that `sort(e_diag, rev=true)[invsort] == e_diag`
            E[1+offset:n_energies+offset, p] .= E[invsort, p]
        end
    end
end

"""
Calculate energy bands of the Floquet Hamiltonian (S20) with boundaries sweeping over the adiabatic `phases` φₓ. 
The operation of this function follows that of [`compute_floquet_bands`](@ref).
Parameter `n` is the number of cells in the lattice; the eigenfunctions will be calculated in the basis of functions sin(𝑗𝑥/𝑛) / √(𝑛π/2).
"""
function compute_floquet_bands_with_boundary(; n::Integer, n_min::Integer, n_max::Integer, phases::AbstractVector{<:Real}, s::Integer, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
    X(j′, j) = 16n*j*j′ / (π*((j-j′)^2-(2n)^2)*((j+j′)^2-(2n)^2))
    
    gs1 = 2n - 1 # number of levels in the first band of spatial Hamiltonian (group size 1)
    gs2 = 2n + 1 # number of levels in the second band of spatial Hamiltonian (group size 2)
    # convert `n_min` and `n_max` to actual level numbers
    n_min = (n_min-1) ÷ 2 * 4n + (isodd(n_min) ? 1 : gs1 + 1)
    n_max = (n_max-1) ÷ 2 * 4n + (isodd(n_max) ? gs1 : gs1 + gs2)
    if iseven(n_min) # swap `gs1` and `gs2` so that they correspond to actual group sizes
        gs1, gs2 = gs2, gs1
    end

    n_j = 2n_max # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian
    h = zeros(n_j, n_j)

    Δn = n_max - n_min + 1
    ν = Vector{Int}(undef, Δn)
    # FIll `ν`: [1 (`gs1` times), 2 (`gs2` times), 3 (`gs1` times), 4 (`gs2` times), ...]
    number = 1
    g = gs1 + gs2
    for i in 0:Δn÷g-1
        ν[g*i+1:g*i+gs1] .= number
        number += 1
        ν[g*i+gs1+1:g*i+g] .= number
        number += 1
    end
    ν[Δn - Δn%g + 1:end] .= number

    pattern = [fill(gs1, gs1); fill(gs2, gs2)]
    G = repeat(pattern, Δn÷g) # a pattern which e.g. for `n == 2` looks like [3, 3, 3, 5, 5, 5, 5, 3, 3, 3, 5, 5, 5, 5, ...]
    Δn % g != 0 && append!(G, fill(gs1, gs1))
    
    ϵ = Matrix{Float64}(undef, Δn, length(phases)) # eigenvalues of ℎ (the unperturbed Hamiltonian)
    c = Matrix{Float64}(undef, n_j, Δn) # eigenvectors of ℎ
    
    E = Matrix{Float64}(undef, Δn, length(phases)) # eigenvalues of 𝐻 (Floquet quasi-energies)
    H_dim = Δn # dimension of the constructed 𝐻 matrix
    # number of non-zero elements in 𝐻:
    n_H_nonzeros = H_dim + 2*( # diagonal plus two times upper off-diagonal terms:
                   (H_dim ÷ g - 1) * (gs1^2 + gs2^2) + # number of long  lattice blocks of size `g` × `g`, each having `(gs1^2 + gs2^2)` elements
                   (H_dim ÷ g - 2) * (gs1^2 + gs2^2) + # number of short lattice blocks of size `g` × `g`, each having `(gs1^2 + gs2^2)` elements
                   2(H_dim % g != 0 ? gs1^2 : 0) ) # if `H_dim % g != 0`, then one more block of size `gs1` is present, both for short and long lattice
   
    H_rows = zeros(Int, n_H_nonzeros)
    H_cols = zeros(Int, n_H_nonzeros)
    H_vals = zeros(ComplexF64, n_H_nonzeros)
    
    for (z, ϕ) in enumerate(phases)
        if pumptype != :time || z == 1 # If pupming is not time-only, ℎ has to be diagonalised on each iteration. If it's time-only, then we diagonalise only once, at `z == 1`.
            for j in 1:n_j
                for j′ in j:n_j
                    val = 0.0
                    if isodd(j′ + j)
                        val += Vₗ/2 * X(j′, j) * sin(2ϕ)
                    else
                        # check diagonals "\"
                        if j′ == j
                            val += (gₗ + Vₗ)/2 + (j / n)^2
                        elseif j′ == j - 2n || j′ == j + 2n
                            val += Vₗ/2 * cos(2ϕ) / 2
                        elseif j′ == j - 4n || j′ == j + 4n
                            val += gₗ/2 / 2
                        end
                        # check anti-diagonals "/"
                        if j′ == -j - 2n || j′ == -j + 2n
                            val += -Vₗ/2 * cos(2ϕ) / 2
                        elseif j′ == -j - 4n || j′ == -j + 4n
                            val += -gₗ/2 / 2
                        end
                    end
                    h[j′, j] = h[j, j′] = val # push the element to the conjugate positions
                end
            end
            f = eigen(h)
            # save only energies and states for levels from `n_min` to `n_max`
            ϵ[:, z] = f.values[n_min:n_max]
            c .= f.vectors[:, n_min:n_max]
            if pumptype == :time
                for p in 2:length(phases) # copy the calculated first column of `ϵ` to all other columns for consistency
                    ϵ[:, p] = ϵ[:, 1]
                end
            end
        end

        # Construct 𝐻
        p = 1 # a counter for placing elements to the vectors `H_*`
        for m in 1:H_dim
            # place the diagonal element (S25)
            H_rows[p] = H_cols[p] = m
            H_vals[p] = ϵ[m, z] - ν[m]*ω/s
            p += 1

            # place the elements of the long lattice (S26)
            for i in 1:G[m]
                # skip `s` groups of `g`, then some more groups depending on `m`, then skip `G[1]` cells
                m′ = g*(s÷2) + g*((ν[m]-1)÷2) + iseven(ν[m])*G[1] + i
                m′ > H_dim && break
                H_rows[p] = m′
                H_cols[p] = m
                if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                    j_sum = sum( (c[j+4n, m′]/4 + c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 4n+1:n_j-4n ) + 
                            sum( (c[j+4n, m′]/4 - c[-j+4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 1:4n-1 ) +
                            (c[4n+4n, m′]/4 + c[4n, m′]/2) * c[4n, m] + # iteration `j = 4n`
                            sum( (c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = n_j-4n+1:n_j )
                    H_vals[p] = (pumptype == :space ? λₗ/2 * j_sum : λₗ/2 * j_sum * cis(-2ϕ)) # a check for space or space-time pumping
                elseif pumptype == :time 
                    H_vals[p] *= cis(-2(phases[2]-phases[1]))
                end
                p += 1
                # place the conjugate element
                H_rows[p] = m
                H_cols[p] = m′
                H_vals[p] = H_vals[p-1]'
                p += 1
            end
            
            # place the elements of the short lattice (S29)
            for i in 1:G[m]
                m′ = g*s + g*((ν[m]-1)÷2) + iseven(ν[m])*G[1] + i
                m′ > H_dim && break
                H_rows[p] = m′
                H_cols[p] = m
                if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                    j_sum = sum( (-c[j+4n, m′]/4 - c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 4n+1:n_j-4n ) + 
                            sum( (-c[j+4n, m′]/4 + c[-j+4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 1:4n-1) +
                            (-c[4n+4n, m′]/4 + c[4n, m′]/2) * c[4n, m] + # iteration `j = 4n`
                            sum( (-c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = n_j-4n+1:n_j)
                    H_vals[p] = λₛ/2 * j_sum
                end
                p += 1
                # place the conjugate element
                H_rows[p] = m
                H_cols[p] = m′
                H_vals[p] = H_vals[p-1]'
                p += 1
            end
        end
        H = sparse(H_rows, H_cols, H_vals)
        vals, _, info = eigsolve(H, Δn, :LR; krylovdim=H_dim)
        if info.converged < Δn
            @warn "Only $(info.converged) eigenvalues out of $(Δn) converged when diagonalising 𝐻ₖ. "*
                  "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
        end
        E[:, z] .= vals[1:Δn]
    end
    return ϵ, E
end

"""
Calculate energy bands of the Floquet Hamiltonian (S20) with boundaries sweeping over the adiabatic `phases` φₓ. 
The operation of this function follows that of [`compute_floquet_bands`](@ref).
Parameter `n` is the number of cells in the lattice; the eigenfunctions will be calculated in the basis of functions sin(𝑗𝑥/𝑛) / √(𝑛π/2).
"""
function compute_floquet_bands_states(; n::Integer, n_min::Integer, n_max::Integer, phases::AbstractVector{<:Real}, s::Integer, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
    X(j′, j) = 16n*j*j′ / (π*((j-j′)^2-(2n)^2)*((j+j′)^2-(2n)^2))
    
    gs1 = 2n - 1 # number of levels in the first band of spatial Hamiltonian (group size 1)
    gs2 = 2n + 1 # number of levels in the second band of spatial Hamiltonian (group size 2)
    # convert `n_min` and `n_max` to actual level numbers
    n_min = (n_min-1) ÷ 2 * 4n + (isodd(n_min) ? 1 : gs1 + 1)
    n_max = (n_max-1) ÷ 2 * 4n + (isodd(n_max) ? gs1 : gs1 + gs2)
    if iseven(n_min) # swap `gs1` and `gs2` so that they correspond to actual group sizes
        gs1, gs2 = gs2, gs1
    end

    n_j = 2n_max # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian
    h = zeros(n_j, n_j)

    Δn = n_max - n_min + 1
    ν = Vector{Int}(undef, Δn)
    # FIll `ν`: [1 (`gs1` times), 2 (`gs2` times), 3 (`gs1` times), 4 (`gs2` times), ...]
    number = 1
    g = gs1 + gs2
    for i in 0:Δn÷g-1
        ν[g*i+1:g*i+gs1] .= number
        number += 1
        ν[g*i+gs1+1:g*i+g] .= number
        number += 1
    end
    ν[Δn - Δn%g + 1:end] .= number

    pattern = [fill(gs1, gs1); fill(gs2, gs2)]
    G = repeat(pattern, Δn÷g) # a pattern which e.g. for `n == 2` looks like [3, 3, 3, 5, 5, 5, 5, 3, 3, 3, 5, 5, 5, 5, ...]
    Δn % g != 0 && append!(G, fill(gs1, gs1))
    
    H_dim = Δn # dimension of the constructed 𝐻 matrix
    # number of non-zero elements in 𝐻:
    n_H_nonzeros = H_dim + 2*( # diagonal plus two times upper off-diagonal terms:
                   (H_dim ÷ g - 1) * (gs1^2 + gs2^2) + # number of long  lattice blocks of size `g` × `g`, each having `(gs1^2 + gs2^2)` elements
                   (H_dim ÷ g - 2) * (gs1^2 + gs2^2) + # number of short lattice blocks of size `g` × `g`, each having `(gs1^2 + gs2^2)` elements
                   2(H_dim % g != 0 ? gs1^2 : 0) ) # if `H_dim % g != 0`, then one more block of size `gs1` is present, both for short and long lattice
   
    H_rows = zeros(Int, n_H_nonzeros)
    H_cols = zeros(Int, n_H_nonzeros)
    H_vals = zeros(ComplexF64, n_H_nonzeros)
    
    ϕ = phases[1]

    for j in 1:n_j
        for j′ in j:n_j
            val = 0.0
            if isodd(j′ + j)
                val += Vₗ/2 * X(j′, j) * sin(2ϕ)
            else
                # check diagonals "\"
                if j′ == j
                    val += (gₗ + Vₗ)/2 + (j / n)^2
                elseif j′ == j - 2n || j′ == j + 2n
                    val += Vₗ/2 * cos(2ϕ) / 2
                elseif j′ == j - 4n || j′ == j + 4n
                    val += gₗ/2 / 2
                end
                # check anti-diagonals "/"
                if j′ == -j - 2n || j′ == -j + 2n
                    val += -Vₗ/2 * cos(2ϕ) / 2
                elseif j′ == -j - 4n || j′ == -j + 4n
                    val += -gₗ/2 / 2
                end
            end
            h[j′, j] = h[j, j′] = val # push the element to the conjugate positions
        end
    end
    f = eigen(h)
    # save only energies and states for levels from `n_min` to `n_max`
    ϵ = f.values[n_min:n_max]
    c = f.vectors[:, n_min:n_max]

    # Construct 𝐻
    p = 1 # a counter for placing elements to the vectors `H_*`
    for m in 1:H_dim
        # place the diagonal element (S25)
        H_rows[p] = H_cols[p] = m
        H_vals[p] = ϵ[m] - ν[m]*ω/s
        p += 1

        # place the elements of the long lattice (S26)
        for i in 1:G[m]
            # skip `s` groups of `g`, then some more groups depending on `m`, then skip `G[1]` cells
            m′ = g*(s÷2) + g*((ν[m]-1)÷2) + iseven(ν[m])*G[1] + i
            m′ > H_dim && break
            H_rows[p] = m′
            H_cols[p] = m
                j_sum = sum( (c[j+4n, m′]/4 + c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 4n+1:n_j-4n ) + 
                        sum( (c[j+4n, m′]/4 - c[-j+4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 1:4n-1 ) +
                        (c[4n+4n, m′]/4 + c[4n, m′]/2) * c[4n, m] + # iteration `j = 4n`
                        sum( (c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = n_j-4n+1:n_j )
                H_vals[p] = (pumptype == :space ? λₗ/2 * j_sum : λₗ/2 * j_sum * cis(-2ϕ)) # a check for space or space-time pumping

            p += 1
            # place the conjugate element
            H_rows[p] = m
            H_cols[p] = m′
            H_vals[p] = H_vals[p-1]'
            p += 1
        end
        
        # place the elements of the short lattice (S29)
        for i in 1:G[m]
            m′ = g*s + g*((ν[m]-1)÷2) + iseven(ν[m])*G[1] + i
            m′ > H_dim && break
            H_rows[p] = m′
            H_cols[p] = m
                j_sum = sum( (-c[j+4n, m′]/4 - c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 4n+1:n_j-4n ) + 
                        sum( (-c[j+4n, m′]/4 + c[-j+4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = 1:4n-1) +
                        (-c[4n+4n, m′]/4 + c[4n, m′]/2) * c[4n, m] + # iteration `j = 4n`
                        sum( (-c[j-4n, m′]/4 + c[j, m′]/2) * c[j, m] for j = n_j-4n+1:n_j)
                H_vals[p] = λₛ/2 * j_sum
            p += 1
            # place the conjugate element
            H_rows[p] = m
            H_cols[p] = m′
            H_vals[p] = H_vals[p-1]'
            p += 1
        end
    end
    H = sparse(H_rows, H_cols, H_vals)
    E, b, info = eigsolve(H, Δn, :LR; krylovdim=H_dim)
    if info.converged < Δn
        @warn "Only $(info.converged) eigenvalues out of $(Δn) converged when diagonalising 𝐻ₖ. "*
                "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
    end

    return ϵ, E, c, b
end

"""
Permute Floquet energy levels calculated with open boundary conditions contained in `E` so that they are stored in the same order as the eigenenergies `e` of
the spatial Hamiltonian.
The operation of this function follows that of [`permute_floquet_bands`](@ref).
"""
function permute_floquet_bands_with_boundary!(E::AbstractMatrix{<:Float64}, e::AbstractMatrix{<:Float64}; n_cells::Integer, n_min::Integer, ω::Real, s::Integer)
    n_energies, n_phases = size(e)

    gs1 = 2n_cells - 1 # number of levels in the first band of spatial Hamiltonian (group size 1)
    gs2 = 2n_cells + 1 # number of levels in the second band of spatial Hamiltonian (group size 2)
    if iseven(n_min) # swap `gs1` and `gs2` so that they correspond to actual group sizes
        gs1, gs2 = gs2, gs1
    end

    ν = Vector{Int}(undef, n_energies)
    # FIll `ν`: [1 (`gs1` times), 2 (`gs2` times), 3 (`gs1` times), 4 (`gs2` times), ...]
    number = 1;
    g = gs1 + gs2
    for i in 0:n_energies÷g-1
        ν[g*i+1:g*i+gs1] .= number
        number += 1
        ν[g*i+gs1+1:g*i+g] .= number
        number += 1
    end
    ν[n_energies - n_energies%g + 1:end] .= number
    ν .*= ω/s
    
    for p in 1:n_phases
        e_diag = [e[m, p] - ν[m] for m in 1:n_energies] # Floquet energies at zero perturbation
        invsort = sortperm(sortperm(e_diag, rev=true))  # inverse permutation, such that `sort(e_diag, rev=true)[invsort] == e_diag`
        E[1:n_energies, p] .= E[invsort, p]
    end
end

