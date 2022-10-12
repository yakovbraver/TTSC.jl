module Bandsolvers

using LinearAlgebra: eigen, eigvals, schur, ⋅, diagm, diagind

"A type for storing the Wannier functions."
mutable struct Wanniers
    minlevel::Int # number of the first energy level of ℎ to use for constructing wanniers (this is used in the unperturbed case)
    targetlevels::Vector{Int} # numbers of quasienergy levels to use for constructing wanniers (this is used in the Floquet case)
    n_lo::Vector{Int} # number of levels in the lower subband at each phase; in non-periodic case this depends on the edge state branch position 
    E::Matrix{Float64} # `E[j, i]` = mean energy of `j`th wannier at `i`th phase; use nested arrays because in non-periodic case
    pos::Matrix{Float64} # `pos[j, i]` = position eigenvalue of `j`th wannier at `i`th phase
    d::Array{ComplexF64, 3} # `d[:, :, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details
end

"Default-construct an empty `Wanniers` object."
function Wanniers()
    Wanniers(0, Int[], Int[], Float64[;;], Float64[;;], ComplexF64[;;;])
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
    w::Wanniers
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

    UnperturbedHamiltonian(Int(n_cells), Float64(M), (l === nothing ? 1 : l), Float64(gₗ), Float64(Vₗ), isperiodic,
                           collect(Float64, phases), maxlevel, bandsizes, E, c, Wanniers())
end

"Diagonalise the unperturbed Hamiltonian `uh` at each phase."
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

"Calculate Wannier vectors for the unperturbed Hamiltonian using the energy eigenstates in the band number `targetband`."
function compute_wanniers!(uh::UnperturbedHamiltonian; targetband::Integer)
    N = uh.N

    minlevel = (targetband-1) * 2N + 1
    uh.w.minlevel = minlevel # save this because it's needed in `make_wavefunction` when constructing coordinate space Wannier functions

    if uh.isperiodic
        E = Matrix{Float64}(undef, 2N, length(uh.phases))
        pos = Matrix{Float64}(undef, 2N, length(uh.phases))
        # `d` fill format: `d[1:N, 1:N, i]` = eigenvectors of the lower subband, `d[1:N, N+1:2N, i]` = eigenvectors of the higher subband
        d = Array{ComplexF64, 3}(undef, N, 2N, length(uh.phases))
        uh.w = Wanniers(minlevel, Int[], fill(N, length(uh.phases)), E, pos, d)

        X = Matrix{ComplexF64}(undef, N, N) # position operator
        
        for i in eachindex(uh.phases)
            for o in (0, N)
                window = 1+o:N+o
                for n in 1:N
                    for n′ in 1:N
                        X[n′, n] = sum(uh.c[j+1, minlevel+o+n′-1, i]' * uh.c[j, minlevel+o+n-1, i] for j = 1:size(uh.c, 1)-1)
                    end
                end
                # `eigen` does not guarantee orthogonality of eigenvectors in case of degeneracies for `X` unitary, so use `schur`
                # (although a degeneracy of coordinates eigenvalues is unlikely here)
                _, d[:, window, i], pos_complex = schur(X)
                pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
                sp = sortperm(pos_real)                 # sort the eigenvalues
                pos[window, i] = pos_real[sp]
                @views Base.permutecols!!(d[:, window, i], sp)    # sort the eigenvectors in the same way
                E[window, i] = [abs2.(dˣ) ⋅ uh.E[minlevel+o:minlevel+N+o-1, i] for dˣ in eachcol(d[:, window, i])]
            end
        end
    else
        n_w = isodd(targetband) ? 2N-1 : 2N+1 # total number of wanniers to construct; this is the number of levels in the target band
        E = Matrix{Float64}(undef, n_w, length(uh.phases))
        pos = Matrix{Float64}(undef, n_w, length(uh.phases))
        # `d` fill format: `d[1:n_lo[i], 1:n_lo[i], i]` = eigenvectors of the lower subband,
        #                  `d[1:n_w-n_lo[i], n_lo[i]+1:n_w, i]` = eigenvectors of the higher subband
        d = Array{ComplexF64, 3}(undef, n_w÷2+1, n_w, length(uh.phases))
        uh.w = Wanniers(minlevel, Int[], Vector{Int}(undef, length(uh.phases)), E, pos, d)

        X_less = zeros(n_w÷2, n_w÷2) # position operator for a subband which does not contain the edge state branch
        X_more = zeros(n_w÷2+1, n_w÷2+1) # position operator for a subband which contains the edge state branch
        
        n_j = size(uh.c, 1)

        for i in eachindex(uh.phases)
            up = uh.E[minlevel+n_w÷2, i] > (uh.E[minlevel+n_w÷2-1, i] + uh.E[minlevel+n_w÷2+1, i])/2 # true if the edge state branch is above the mean value
            
            # Lower band
            X = up ? X_less : X_more # bind the position operator to a matrix of the appropriate size
            n_lo = n_w÷2 + !up # number of levels in the lower subband
            uh.w.n_lo[i] = n_lo
            for n in 1:n_lo
                for n′ in n:n_lo
                    X[n′, n] = X[n, n′] = N*sum(uh.c[j, minlevel+n-1, i] * (π/2 * uh.c[j, minlevel+n′-1, i] - 8/π * sum(uh.c[j′, minlevel+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                            for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
                end
            end
            uh.w.pos[1:n_lo, i], uh.w.d[1:n_lo, 1:n_lo, i] = eigen(X)
            uh.w.E[1:n_lo, i] = [dˣ.^2 ⋅ uh.E[minlevel:minlevel+n_lo-1, i] for dˣ in eachcol(uh.w.d[1:n_lo, 1:n_lo, i])]

            # Higher band
            X = up ? X_more : X_less
            n_hi = n_w - n_lo
            for n in 1:n_hi
                for n′ in n:n_hi
                    X[n′, n] = X[n, n′] = N*sum(uh.c[j, minlevel+n_lo+n-1, i] * (π/2 * uh.c[j, minlevel+n_lo+n′-1, i] - 8/π * sum(uh.c[j′, minlevel+n_lo+n′-1, i]*j*j′/(j^2-j′^2)^2
                                                                                 for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
                end
            end
            uh.w.pos[n_lo+1:n_w, i], uh.w.d[1:n_hi, n_lo+1:n_w, i] = eigen(X)
            uh.w.E[n_lo+1:n_w, i] = [dˣ.^2 ⋅ uh.E[minlevel+n_lo:minlevel+n_w-1, i] for dˣ in eachcol(uh.w.d[1:n_hi, n_lo+1:n_w, i])]
        end
    end
end

"""
Construct energy eigenfunctions at coordinates in `x` for each eigenstate number in `whichstates` at each phase number in `whichphases`.
Return `ψ`, where `ψ[:, j, i]` = `j`th eigenfunction at `i`th phase.
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
Construct Wannier functions at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `uh` are constructed.
In the process, energy eigenfunctions are also constructed.
Return `ψ, w`, where `ψ[:, j, i]` = `j`th eigenfunction at `i`th phase, and `w[:, j, i]` = `j`th Wannier function at `i`th phase.
"""
function make_wannierfunctions(uh::UnperturbedHamiltonian, x::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
    n_w = size(uh.w.E, 1)
    w = Array{ComplexF64, 3}(undef, length(x), n_w, length(uh.phases))
    ψ = make_eigenfunctions(uh, x, whichphases, uh.w.minlevel:uh.w.minlevel+n_w) # construct energy eigenfunctions
    for i in eachindex(whichphases)
        for j in 1:uh.w.n_lo[i]
            w[:, j, i] = sum(uh.w.d[k, j, i] * ψ[:, uh.w.minlevel+k-1, i] for k = 1:uh.w.n_lo[i])
        end
        for j in uh.w.n_lo[i]+1:n_w
            w[:, j, i] = sum(uh.w.d[k, j, i] * ψ[:, uh.w.minlevel+uh.w.n_lo[i]+k-1, i] for k = 1:n_w-uh.w.n_lo[i])
        end
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
A type representing the Floquet Hamiltonian
    ℋ = ℎ - i∂ₜ + λₛsin²(2𝑥)cos(2𝜔𝑡) + λₗcos²(2𝑥)cos(𝜔𝑡 + 𝜑ₜ),
where ℎ is the unperturbed Hamiltonian represented by [`UnperturbedHamiltonian`](@ref).
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

"""
Construct a `FloquetHamiltonian` object. `minband` is the first energy band of `uh` to use when constructing the Floquet Hamiltonian matrix.
Type of pumping is controlled via `pumptype`: `:time` for temporal, `:space` for spatial, or anything else for simultaneous space-time pumping.
"""
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

"Diagonalise the Floquet Hamiltonian `fh` at each phase."
function diagonalise!(fh::FloquetHamiltonian)
    # make views
    (;N, phases, E, c) = fh.uh
    n_j = size(c, 1)
    (;s, ω, λₛ, λₗ, pumptype) = fh

    n_levels = fh.uh.maxlevel - fh.minlevel + 1 # number of levels of spatial Hamiltonian to use for constructing Floquet Hamiltonian
    ν(m) = ceil(Int, m/2N)

    H = zeros(ComplexF64, n_levels, n_levels) # Floquet Hamiltonian matrix

    for (i, ϕ) in enumerate(fh.uh.phases)
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

"""
Permute Floquet quasienergy levels contained in `fh.E` so that they are stored in the same order as the eigenenergies of ℎ stored in `fh.uh.E`.
Repeat this for every phase (i.e. column of `fh.E`).
To perfrorm the sorting, we first calculate `fh.uh.E - ν(m)` which is the diagonal of ℋ. If there is no perturbation, then these
are the Floquet quasienergies. Then, we sort them in descending order (as if we diagonalised the Hamiltonian) and find the permutation
that would undo this sorting. This permutation is applied to a copy of `fh.E`.
The procedure yields fully correct results only if `E` has been calculated at zero perturbation. The perturbation may additionally change
the order of levels, and there is no simple way to disentangle the order. The permutation is still useful in that case, but the results 
should not be taken too literally.
"""
function order_floquet_levels(fh::FloquetHamiltonian)
    E = similar(fh.E)
    ν(m) = ceil(Int, m/2fh.uh.N)
    n_levels = size(fh.E, 1) # number of Floquet levels; equivalently, number of levels of ℎ used for constructing ℋ
    for i in eachindex(fh.uh.phases)
        E_diag = [fh.uh.E[m, i] - ν(m) * fh.ω/fh.s for m in 1:n_levels] # Floquet energies at zero perturbation
        invsort = sortperm(sortperm(E_diag, rev=true))  # inverse permutation, such that `sort(E_diag, rev=true)[invsort] == E_diag`
        E[:, i] .= fh.E[invsort, i]
    end
    return E
end

# function compute_wanniers!(fh::FloquetHamiltonian; targetband::Integer)
#     # n_target_min = (n_target-1) * 4N + 1

#     # n_j = n_max * 2N # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian

#     # h = zeros(ComplexF64, 2n_j + 1, 2n_j + 1)
#     # h = diagm(0 => ComplexF64[(2j/N)^2 + (gₗ + Vₗ)/2 for j = -n_j:n_j])
#     # h[diagind(h, -2N)] .= h[diagind(h, 2N)] .= gₗ/4

#     # n_w = mix_time_cells ? s*N : N # number of Wannier functions to construct
#     # pos_lo = Matrix{Float64}(undef, n_w, length(phases)) # position eigenvalues (Wannier centres) of the lower spatial levels
#     # pos_hi = Matrix{Float64}(undef, n_w, length(phases)) # position eigenvalues (Wannier centres) of the higher spatial levels
#     # ε_lo = Matrix{Float64}(undef, n_w, length(phases)) # energies of Wannier states of the lower spatial levels
#     # ε_hi = Matrix{Float64}(undef, n_w, length(phases)) # energies of Wannier states of the higher spatial levels
#     # wf_lo = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the lower spatial levels
#     # wf_hi = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the higher spatial levels
#     # window_lo = Int[] # window of Floquet eigenstates which will be used to construct Wannier functions out of the lower spatial levels
#     # window_hi = Int[] # window of Floquet eigenstates which will be used to construct Wannier functions out of the higher spatial levels
#     # for i in 0:n_w÷N - 1
#     #     append!(window_hi, n_target_min+2i*N:n_target_min+2i*N + N - 1)
#     #     append!(window_lo, n_target_min+(2i+1)*N:n_target_min+(2i+1)*N + N - 1)
#     # end

#     # u_lo = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the lower spatial levels
#     # u_hi = Array{Float64,4}(undef, length(coords), n_w, length(ωts), length(phases)) # Wannier states of the higher spatial levels

#     # make views
#     (;N, phases, E, c) = fh.uh
#     n_j = size(c, 1)
#     (;s, ω, λₛ, λₗ, pumptype) = fh

#     n_levels = fh.uh.maxlevel - fh.minlevel + 1 # number of levels of spatial Hamiltonian to use for constructing Floquet Hamiltonian
#     ν(m) = ceil(Int, m/2N)

#     # ψ = Matrix{ComplexF64}(undef, length(coords), n_levels) # eigenvectors of ℎ in 𝑥-representation
#     # cc = Matrix{ComplexF64}(undef, n_levels, n_levels) # matrix of products of `c`'s that will be needed multiple times
#     # ccc = Matrix{ComplexF64}(undef, n_levels, n_levels) # products `cc`'s and cis
    
#     # x = Matrix{ComplexF64}(undef, n_w, n_w) # position operator

#     H = zeros(ComplexF64, n_levels, n_levels) # Floquet Hamiltonian matrix

#     for (i, ϕ) in enumerate(fh.uh.phases)
#             # construct coordinate representation of eigenfunctions and compute products of `c`'s that will be needed multiple times
#             for m in 1:n_levels
#                 ψ[:, m] = make_exp_state(coords, c[:, m]; n=N)
#                 for m′ in 1:n_levels
#                     cc[m′, m] = sum(c[j+1, m′]' * c[j, m] for j = 1:2n_j)
#                 end
#             end

#         t = (pumptype == :space ? π/5 : π/5 - z/length(phases)*π/2) # time moment at which to diagonalise the coordinate operator
#         for m in 1:n_levels, m′ in 1:n_levels
#             ccc[m′, m] = cc[m′, m] * cis((ν(m′)-ν(m))*t)
#         end

#         # Higher band
#         # the loop below runs faster if we make a copy rather than a view of `f.vectors`; 
#         # both approaches are ~6 times faster compared to iterating directly over `f.vectors`
#         b .= f.vectors[:, window_hi]
#         for n in 1:n_w, n′ in 1:n_w
#             x[n′, n] = sum(b[m, n] * sum(b[m′, n′]' * ccc[m′, m] for m′ in 1:n_levels) for m in 1:n_levels)
#         end
#         _, d, pos_complex = schur(x)
#         pos_real = (angle.(pos_complex) .+ π) / 2π * N*π # take angle, convert from (-π, π) to (0, 2π), and map to the interval (0, Nπ)
#         sp = sortperm(pos_real)
#         pos_hi[:, z] = pos_real[sp]   # sort positions
#         Base.permutecols!!(d, sp)     # sort eigenvectors in the same way
#         ε_hi[:, z] = [abs2.(dˣ) ⋅ E[window_hi, z] for dˣ in eachcol(d)]
#         for (t, ωt) in enumerate(ωts)
#             for X in 1:n_w
#                 wf_hi[:, X, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * sum(d[l, X] * b[m, l] for l = 1:n_w) for m in 1:n_levels))
#             end
#         end

#         # Lower band
#         b .= f.vectors[:, window_lo]
#         for n in 1:n_w, n′ in 1:n_w
#             x[n′, n] = sum(b[m, n] * sum(b[m′, n′]' * ccc[m′, m] for m′ in 1:n_levels) for m in 1:n_levels)
#         end
#         _, d, pos_complex = schur(x)
#         pos_real = (angle.(pos_complex) .+ π) / 2π * N*π # take angle, convert from (-π, π) to (0, 2π), and map to the interval (0, Nπ)
#         sp = sortperm(pos_real)
#         pos_lo[:, z] = pos_real[sp]   # sort positions
#         Base.permutecols!!(d, sp)     # sort eigenvectors in the same way
#         ε_lo[:, z] = [abs2.(dˣ) ⋅ E[window_lo, z] for dˣ in eachcol(d)]
#         for (t, ωt) in enumerate(ωts)
#             for X in 1:n_w
#                 wf_lo[:, X, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * sum(d[l, X] * b[m, l] for l = 1:n_w) for m in 1:n_levels))
#             end
#         end
#     end
# end

"""
Calculate Wannier vectors for the FloquetHamiltonian Hamiltonian `fh` using the quasienergy levels `targetlevels_lo` and `targetlevels_hi`.
"""
function compute_wanniers!(fh::FloquetHamiltonian; targetlevels_lo::AbstractVector{<:Real}, targetlevels_hi::AbstractVector{<:Real})
    (;N, phases, c) = fh.uh
    (;pos_lo, pos_hi, E_lo, E_hi, d_lo, d_hi) = fh.uh.w
    ν(m) = ceil(Int, m/2N)

    fh.uh.w.targetlevels_lo = targetlevels_lo # save this because it's needed in `make_wavefunction` when constructing coordinate space Wannier functions
    fh.uh.w.targetlevels_hi = targetlevels_hi

    X = Matrix{ComplexF64}(undef, length(targetlevels_lo), length(targetlevels_lo)) # position operator
    
    n_levels = fh.uh.maxlevel - fh.minlevel + 1 # number of levels of spatial Hamiltonian to use for constructing Floquet Hamiltonian
    ∑cc = Matrix{ComplexF64}(undef, n_levels, n_levels)
    cis∑cc = Matrix{ComplexF64}(undef, n_levels, n_levels)
    
    for i in eachindex(phases)
        for m in fh.minlevel:fh.minlevel+n_levels-1
            for m′ in fh.minlevel:fh.minlevel+n_levels-1
                ∑cc[m′-fh.minlevel+1, m-fh.minlevel+1] = sum(c[j+1, m′, i]' * c[j, m, i] for j = 1:size(c, 1)-1)
            end
        end

        t = (fh.pumptype == :space ? π/5 : π/5 - i/length(phases)*π/2) # time moment at which to diagonalise the coordinate operator
        if fh.pumptype == :space || i == 1 # if pumping is space-only, `cis∑cc` can be calculated only once, at `i == 1`, because `t` does not change
            for m in 1:n_levels
                for m′ in 1:n_levels
                    cis∑cc[m′, m] = ∑cc[m′, m] * cis((ν(m′+fh.minlevel-1) - ν(m+fh.minlevel-1)) * t)
                end
            end
        end

        # Lower band
        for (in, n) in enumerate(targetlevels_lo)
            for (in′, n′) in enumerate(targetlevels_lo)
                X[in′, in] = sum(fh.b[m, n, i] * sum(fh.b[m′, n′, i]' * cis∑cc[m′, m] for m′ in 1:n_levels) for m in 1:n_levels)
            end
        end
        # `eigen` does not guarantee orthogonality of eigenvectors in case of degeneracies for `X` unitary, so use `schur`
        # (although a degeneracy of coordinates eigenvalues is unlikely here)
        _, d_lo[i], pos_complex = schur(X)
        pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
        sp = sortperm(pos_real)         # sort the eigenvalues
        pos_lo[i] = pos_real[sp]
        Base.permutecols!!(d_lo[i], sp) # sort the eigenvectors in the same way
        E_lo[i] = [abs2.(dˣ) ⋅ fh.E[targetlevels_lo, i] for dˣ in eachcol(d_lo[i])]

        # Higher band
        for (in, n) in enumerate(targetlevels_hi)
            for (in′, n′) in enumerate(targetlevels_hi)
                X[in′, in] = sum(fh.b[m, n, i] * sum(fh.b[m′, n′, i]' * cis∑cc[m′, m] for m′ in 1:n_levels) for m in 1:n_levels)
            end
        end
        _, d_hi[i], pos_complex = schur(X)
        pos_real = @. mod2pi(angle(pos_complex)) / 2 * N # `mod2pi` converts the angle from [-π, π) to [0, 2π)
        sp = sortperm(pos_real)
        pos_hi[i] = pos_real[sp]
        Base.permutecols!!(d_hi[i], sp)
        E_hi[i] = [abs2.(dˣ) ⋅ fh.E[targetlevels_hi, i] for dˣ in eachcol(d_hi[i])]
    end
end

# """
# Construct Wannier functions at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `fh` are constructed.
# Return `(w_lo, w_hi)`, where `w_**[i][j][ix, it]` = `j`th Wannier function at `i`th phase at `ix`th coordinate at `it`th time moment.
# """
# function make_wannierfunctions(fh::FloquetHamiltonian, x::AbstractVector{<:Real}, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
#     u = make_eigenfunctions(fh, x, Ωt, whichphases, [fh.uh.w.targetlevels_lo; fh.uh.w.targetlevels_hi])
#     w_lo = [Vector{Matrix{ComplexF64}}() for _ in eachindex(whichphases)]
#     w_hi = [Vector{Matrix{ComplexF64}}() for _ in eachindex(whichphases)]
#     for (i, iϕ) in enumerate(whichphases)
#         n_lo = size(fh.uh.w.d_lo[iϕ], 1)
#         w_lo[i] = [Matrix{ComplexF64}(undef, length(x), length(Ωt)) for _ in 1:n_lo]
#         for j in 1:n_lo
#             w_lo[i][j] = sum(fh.uh.w.d_lo[iϕ][k, j] * u[:, :, k, iϕ] for k = 1:n_lo)
#         end
#         n_hi = size(fh.uh.w.d_hi[iϕ], 1)
#         w_hi[i] = [Matrix{ComplexF64}(undef, length(x), length(Ωt)) for _ in 1:n_hi]
#         for j in 1:n_hi
#             w_hi[i][j] = sum(fh.uh.w.d_lo[iϕ][k, j] * u[:, :, k+n_lo, iϕ] for k = 1:n_hi)
#         end
#     end
#     return w_lo, w_hi
# end

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