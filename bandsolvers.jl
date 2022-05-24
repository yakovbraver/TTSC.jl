import BandedMatrices as BM
using SparseArrays: sparse
using KrylovKit: eigsolve
using LinearAlgebra: eigen, eigvals, schur, ⋅, diagm, diagind
using ProgressMeter

"""
Calculate `n_bands` of energy bands of Hamiltonian (S32) assuming infinite crystal with a quasimomentum 𝑞,
    𝐻 = (𝑝 + 𝑞)²/2𝑀 + 𝜆ₗ𝐴ₗcos(𝑠𝑥 - 𝜒ₗ - φₜ) + 𝜆ₛ𝐴ₛcos(2𝑠𝑥 - 𝜒ₛ)
on 𝑥 ∈ (0; π) sweeping over the adiabatic `phases` φₜ. Boundary conditions are periodic, hence the basis exp(i𝑗𝑥) / √2π is used.
In the returned matrix of bands, columns enumerate the adiabatic phases, while rows enumerate eigenvalues.
Rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
Rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case λₗAₗcos(sϑ+φₜ) leads to 𝑘 = s/2.
"""
function compute_qc_bands(; n_bands::Integer, phases::AbstractVector{<:Real}, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real, χₛ::Real, χₗ::Real)
    n_j = 2n_bands # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
    # Hamiltonian matrix
    H = BM.BandedMatrix{ComplexF64}(undef, (2n_j + 1, 2n_j + 1), (2, 2))
    H[BM.band(-2)] .= λₛAₛ * cis(χₛ)
    H[BM.band(2)]  .= λₛAₛ * cis(-χₛ)
    
    bands = Matrix{Float64}(undef, 2n_bands, length(phases))
    for k in [0, s÷2] # iterate over the centre of BZ and then the boundary
        H[BM.band(0)] .= [(2j + k)^2 / M for j = -n_j:n_j]
        # `a` and `b` control where to place the eigenvalues depedning on `k`; see function docstring
        a = (k > 0)*n_bands + 1 
        b = a+n_bands - 1
        for (i, ϕ) in enumerate(phases)
            H[BM.band(-1)] .= λₗAₗ * cis( χₗ - ϕ)
            H[BM.band(1)]  .= λₗAₗ * cis(-χₗ + ϕ)
            vals, _, _ = eigsolve(H, n_bands, :LR; krylovdim=n_bands+10)
            bands[a:b, i] .= vals[1:n_bands]
        end
    end
    return bands / 2 # restore the omitted factor
end

"""
Calculate `n_levels` of energy levels of Hamiltonian (S32):
    𝐻 = 𝑝²/2𝑀 + 𝜆ₗ𝐴ₗcos(𝑠𝑥 + 𝜒ₗ - φₜ) + 𝜆ₛ𝐴ₛcos(2𝑠𝑥 + 𝜒ₛ)
on 𝑥 ∈ [0; 2π) sweeping over the adiabatic `phases` φₜ. Boundary conditions are periodic, hence the basis exp(i𝑗𝑥) / √2π is used.
In the returned matrix of levels, columns enumerate the adiabatic phases, while rows enumerate eigenvalues.
The eigenvectors are returned as a triple array: `eigvecs[p][n]` holds an eigenvector of `n`th eigenvalue at `p`th phase.
"""
function compute_qc_bands_pbc(; n_levels::Integer, phases::AbstractVector{<:Real}, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real, χₛ::Real, χₗ::Real)
    n_j = 2n_levels # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    # # Hamiltonian matrix
    H = BM.BandedMatrix(BM.Zeros{ComplexF64}(2n_j + 1, 2n_j + 1), (2s, 2s))
    H[BM.band(0)] .= [j^2 / M for j = -n_j:n_j]
    H[BM.band(-2s)] .= λₛAₛ * cis(χₛ)
    H[BM.band(2s)]  .= λₛAₛ * cis(-χₛ)

    levels = Matrix{Float64}(undef, n_levels, length(phases))
    eigvecs = [[Vector{ComplexF64}(undef, 2n_j+1) for _ in 1:n_levels] for _ in 1:length(phases)]
    for (i, ϕ) in enumerate(phases)
        H[BM.band(-s)] .= λₗAₗ * cis( χₗ - ϕ)
        H[BM.band(s)]  .= λₗAₗ * cis(-χₗ + ϕ)
        vals, vecs,  _ = eigsolve(H, n_levels, :LR; krylovdim=2n_levels)
        levels[:, i] = vals[1:n_levels]
        eigvecs[i] .= vecs[1:n_levels]
    end
    return levels / 2, eigvecs
end

"""
Calculate `n_levels` of energy levels of Hamiltonian (S32):
    𝐻 = 𝑝²/2𝑀 + 𝜆ₗ𝐴ₗcos(𝑠𝑥 + 𝜒ₗ - φₜ) + 𝜆ₛ𝐴ₛcos(2𝑠𝑥 + 𝜒ₛ)
on 𝑥 ∈ [0; 2π) sweeping over the adiabatic `phases` φₜ. Boundary conditions are open, hence the basis sin(𝑗𝑥/𝑛) / √(𝑛π/2) is used.
Parameter `n` is the number of cells in the lattice; 𝑗 runs from 0 to `5n_bands`.
Return a tuple (`bands`, `states`): `bands[:, p]` stores eigenenergies at `p`th phase, while `states[p][:, m]` stores `m`th eigenvector at `p`th phase.
Bands and states are sorted in energy-descending order so that for `M` negative, the bands of interest will be the first ones.
`n_bands` is the number of bands of interest, but a larger Hamiltonian matrix is constructed (of size `5n_bands` × `5n_bands`)
so that the bands of interest are calculated correctly. All `5n_bands` energy levels and eigenstates are returned.
"""
function compute_qc_bands_obc(; n_levels::Integer, phases::AbstractVector{<:Real}, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real, χₗ::Real, χₛ::Real)
    X(j′, j, s) = 16s*j*j′ / (π*((j-j′)^2-4s^2)*((j+j′)^2-4s^2))
    
    n_j = 5n_levels # number of indices 𝑗 to use for constructing the Hamiltonian
    H = zeros(n_j, n_j)
    # for storing eigenstates and eigenvectors, see function docstring for format
    bands = Matrix{Float64}(undef, n_j, length(phases))
    states = [Matrix{Float64}(undef, n_j, n_j) for _ in 1:length(phases)]
    for (i, ϕ) in enumerate(phases)
        for j in 1:n_j
            for j′ in j:n_j
                val = 0.0
                if isodd(j′ + j)
                    val += λₗAₗ*X(j′, j, s)*sin(χₗ - ϕ) + λₛAₛ*X(j′, j, 2s)*sin(χₛ)
                else
                    # check diagonals "\"
                    if j′ == j
                        val += j^2 / 8M
                    elseif j′ == j - 2s || j′ == j + 2s
                        val += λₗAₗ * cos(χₗ - ϕ) / 2
                    elseif j′ == j - 4s || j′ == j + 4s
                        val += λₛAₛ * cos(χₛ) / 2
                    end
                    # check anti-diagonals "/"
                    if j′ == -j + 2s
                        val += -λₗAₗ * cos(χₗ - ϕ) / 2
                    elseif j′ == -j + 4s
                        val += -λₛAₛ * cos(χₛ) / 2
                    end
                end
                H[j′, j] = H[j, j′] = val # push the element to the conjugate positions
            end
        end
        bands[:, i], states[i] = eigen(H, sortby=-) # sort in descending order
    end
    return bands, states
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

"""
Diagonalise Floquet Hamiltonian for a periodic system and calculate Wannier centres.
    `n_min` - lowest band number of spatial Hamiltonian to use when constructing Floquet Hamiltonian
    `n_max` - highest band number of spatial Hamiltonian to consider
    `n_target` - number of Floquet band (counting from the highest) to use for Wannier calculations
"""
function compute_floquet_wannier_centres(; N::Integer, n_min::Integer=1, n_target::Integer, n_max::Integer, phases::AbstractVector{<:Real}, s::Integer, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
    n_target_min = (n_target-1) * 4N + 1
    n_target_max = n_target_min + 4N - 1

    n_j = n_max * 2N # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian

    h = zeros(ComplexF64, 2n_j + 1, 2n_j + 1)
    h = diagm(0 => ComplexF64[(2j/N)^2 + (gₗ + Vₗ)/2 for j = -n_j:n_j])
    h[diagind(h, -2N)] .= h[diagind(h, 2N)] .= gₗ/4

    coords = range(0, N*pi, length=40N) # x's for wavefunctions
    ωts = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
    pos_lower = Matrix{Float64}(undef, 2N, length(phases))
    pos_higher = Matrix{Float64}(undef, 2N, length(phases))
    ε_lower = Matrix{Float64}(undef, 2N, length(phases))
    ε_higher = Matrix{Float64}(undef, 2N, length(phases))
    wf_lower = Array{Float64,4}(undef, length(coords), 2N, length(ωts), length(phases))
    wf_higher = Array{Float64,4}(undef, length(coords), 2N, length(ωts), length(phases))

    n_min = (n_min-1) * 2N + 1 # convert `n_min` to actual level number
    n_max = n_max * 2N # convert `n_max` to actual level number
    Δn = n_max - n_min + 1 # number of levels of spatial Hamiltonian to use for constructing Floquet Hamiltonian
    ν(m) = ceil(Int, m/2N)

    ϵ = Matrix{Float64}(undef, Δn, length(phases)) # eigenvalues of ℎ (the unperturbed Hamiltonian)
    c = Matrix{ComplexF64}(undef, 2n_j+1, Δn) # eigenvectors of ℎ in 𝑗-representation
    ψ = Matrix{ComplexF64}(undef, length(coords), Δn) # eigenvectors of ℎ in 𝑥-representation
    cc = Matrix{ComplexF64}(undef, Δn, Δn) # matrix of products of `c`'s that will be needed multiple times
    ccc = Matrix{ComplexF64}(undef, Δn, Δn) # products `cc`'s and cis
    
    x = Matrix{ComplexF64}(undef, 2N, 2N) # position operator

    H_dim = Δn # dimension of the constructed 𝐻 matrix
    H = zeros(ComplexF64, H_dim, H_dim)
    E = Matrix{Float64}(undef, H_dim, length(phases)) # eigenvalues of 𝐻 (Floquet quasi-energies)
    b = Matrix{ComplexF64}(undef, H_dim, 2N) # matrix of eigenvectors of 𝐻

    @showprogress for (z, ϕ) in enumerate(phases)
        if pumptype != :time || z == 1 # If pupming is not time-only, ℎ has to be diagonalised on each iteration. If it's time-only, then we diagonalise only once, at `z == 1`.
            h[diagind(h, -N)] .= Vₗ/4 * cis(2ϕ)
            h[diagind(h, N)]  .= Vₗ/4 * cis(-2ϕ)
            f = eigen(h)
            # save only energies and states for levels from `n_min` to `n_max`
            ϵ[:, z] = f.values[n_min:n_max]
            c .= f.vectors[:, n_min:n_max]
            # construct coordinate representation of eigenfunctions and compute products of `c`'s that will be needed multiple times
            for m in 1:Δn
                ψ[:, m] = make_exp_state(coords, c[:, m]; n=N)
                for m′ in 1:Δn
                    cc[m′, m] = sum(c[j+1, m′]' * c[j, m] for j = 1:2n_j)
                end
            end
            if pumptype == :time
                for p in 2:length(phases) # copy the calculated first column of `ϵ` to all other columns for consistency
                    ϵ[:, p] = ϵ[:, 1]
                end
            end
        end

        # Construct 𝐻
        for m in 1:H_dim
            # place the diagonal element (S25)
            H[m, m] = ϵ[m, z] - ν(m)*ω/s

            # place the elements of the long lattice (S26)
            for i in 1:2N
                m′ = 2N*(s + ν(m) - 1) +  i
                m′ > H_dim && break
                if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                    j_sum = sum( (                c[j, m′]/2 + c[j+2N, m′]/4)' * c[j, m] for j = 1:2N ) +
                            sum( (c[j-2N, m′]/4 + c[j, m′]/2 + c[j+2N, m′]/4)' * c[j, m] for j = 2N+1:(2n_j+1)-2N ) + 
                            sum( (c[j-2N, m′]/4 + c[j, m′]/2                )' * c[j, m] for j = (2n_j+1)-2N+1:(2n_j+1) )
                    H[m′, m] = (pumptype == :space ? λₗ/2 * j_sum : λₗ/2 * j_sum * cis(-2ϕ)) # a check for space or space-time pumping
                elseif pumptype == :time 
                    H[m′, m] *= cis(-2(phases[2]-phases[1]))
                end
                # place the conjugate element
                H[m, m′] = H[m′, m]'
            end
            
            # place the elements of the short lattice (S29)
            for i in 1:2N
                m′ = 2N*(2s + ν(m) - 1) + i
                m′ > H_dim && break
                if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                    j_sum = sum( (                 c[j, m′]/2 - c[j+2N, m′]/4)' * c[j, m] for j = 1:2N ) +
                            sum( (-c[j-2N, m′]/4 + c[j, m′]/2 - c[j+2N, m′]/4)' * c[j, m] for j = 2N+1:(2n_j+1)-2N ) + 
                            sum( (-c[j-2N, m′]/4 + c[j, m′]/2                )' * c[j, m] for j = (2n_j+1)-2N+1:(2n_j+1) )
                    H[m′, m] = λₛ/2 * j_sum
                end
                # place the conjugate element
                H[m, m′] = H[m′, m]'
            end
        end
        f = eigen(H, sortby=x->-real(x))
        E[:, z] .= real.(f.values[1:Δn]) # save all Floquet quasienergies for plotting the spectrum

        ### Wannier centres
        
        t = (pumptype == :space ? pi/5 : pi/5 - z/length(phases)*pi/2) # time moment at which to diagonalise the coordinate operator
        for m in 1:Δn, m′ in 1:Δn
            ccc[m′, m] = cc[m′, m] * cis((ν(m′)-ν(m))*t)
        end

        # Higher band
        # the loop below runs faster if we make a copy rather than a view of `f.vectors`; 
        # both approaches are ~6 times faster compared to iterating directly over `f.vectors`
        window = [n_target_min:n_target_min + N - 1; n_target_min+2N:n_target_min+2N + N - 1]
        b .= f.vectors[:, window]
        for n in 1:2N, n′ in 1:2N
            x[n′, n] = sum(b[m, n] * sum(b[m′, n′]' * ccc[m′, m] for m′ in 1:Δn) for m in 1:Δn)
        end
        _, d, pos_complex = schur(x)
        pos_real = (angle.(pos_complex) .+ π) / 2π * N*π # take angle and convert from (-π, π) to (0, 2π)
        sp = sortperm(pos_real)
        pos_higher[:, z] = pos_real[sp]   # sort positions
        Base.permutecols!!(d, sp)         # sort eigenvectors in the same way
        ε_higher[:, z] = [abs2.(dˣ) ⋅ E[window, z] for dˣ in eachcol(d)]
        for (t, ωt) in enumerate(ωts)
            for X in 1:2N
                wf_higher[:, X, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * sum(d[l, X] * b[m, l] for l = 1:2N) for m in 1:Δn))
            end
        end

        # Lower band
        window = [n_target_min+N:n_target_min+2N-1; n_target_min+3N:n_target_max]
        b .= f.vectors[:, window]
        for n in 1:2N, n′ in 1:2N
            x[n′, n] = sum(b[m, n] * sum(b[m′, n′]' * ccc[m′, m] for m′ in 1:Δn) for m in 1:Δn)
        end
        _, d, pos_complex = schur(x)
        pos_real = (angle.(pos_complex) .+ π) / 2π * N*π # take angle and convert from (-π, π) to (0, 2π)
        sp = sortperm(pos_real)
        pos_lower[:, z] = pos_real[sp]   # sort positions
        Base.permutecols!!(d, sp)         # sort eigenvectors in the same way
        ε_lower[:, z] = [abs2.(dˣ) ⋅ E[window, z] for dˣ in eachcol(d)]
        for (t, ωt) in enumerate(ωts)
            for X in 1:2N
                wf_lower[:, X, t, z] = abs2.(sum(cis(-ν(m)*ωt) * ψ[:, m] * sum(d[l, X] * b[m, l] for l = 1:2N) for m in 1:Δn))
            end
        end
    end
    return ϵ, E, pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher
end

###
"""
Diagonalise spatial Hamiltonian and calculate Wannier centres in non-periodic case
"""
function compute_wannier_centres(; N::Integer, n_min::Integer, n_max::Integer, n_target::Integer, phases::AbstractVector{<:Real}, s::Integer, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real)
    X(j′, j) = 16N*j*j′ / (π*((j-j′)^2-(2N)^2)*((j+j′)^2-(2N)^2))
    
    gs1 = 2N - 1 # number of levels in the first band of spatial Hamiltonian (group size 1)
    gs2 = 2N + 1 # number of levels in the second band of spatial Hamiltonian (group size 2)
    # convert `n_min` and `n_max` to actual level numbers
    n_min = (n_min-1) ÷ 2 * 4N + (isodd(n_min) ? 1 : gs1 + 1)
    n_max = (n_max-1) ÷ 2 * 4N + (isodd(n_max) ? gs1 : gs1 + gs2)
  
    n_target_min = (n_target-1) ÷ 2 * 4N + (isodd(n_target) ? 1 : gs1 + 1)
    n_target_max = n_target_min + (isodd(n_target_min) ? gs1 : gs2) - 1
  

    n_j = 2n_max # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian
    h = zeros(n_j, n_j)

    coords = range(0, N*pi, length=50N) # x's for wavefunctions
    n_w = n_target_max - n_target_min + 1 # numebr of Wannier levels
    pos_lower = [Float64[] for _ in 1:length(phases)]
    pos_higher = [Float64[] for _ in 1:length(phases)]
    ε_lower = [Float64[] for _ in 1:length(phases)]
    ε_higher = [Float64[] for _ in 1:length(phases)]
    wf_lower = [zeros(length(coords), n_w÷2 + 1) for _ in 1:length(phases)]
    wf_higher = [zeros(length(coords), n_w÷2 + 1) for _ in 1:length(phases)]

    x = zeros(n_w÷2 + 1, n_w÷2 + 1)
    
    for (z, ϕ) in enumerate(phases)
        for j in 1:n_j
            for j′ in j:n_j
                val = 0.0
                if isodd(j′ + j)
                    val += Vₗ/2 * X(j′, j) * sin(2ϕ)
                else
                    # check diagonals "\"
                    if j′ == j
                        val += (gₗ + Vₗ)/2 + (j / N)^2
                    elseif j′ == j - 2N || j′ == j + 2N
                        val += Vₗ/2 * cos(2ϕ) / 2
                    elseif j′ == j - 4N || j′ == j + 4N
                        val += gₗ/2 / 2
                    end
                    # check anti-diagonals "/"
                    if j′ == -j - 2N || j′ == -j + 2N
                        val += -Vₗ/2 * cos(2ϕ) / 2
                    elseif j′ == -j - 4N || j′ == -j + 4N
                        val += -gₗ/2 / 2
                    end
                end
                h[j′, j] = h[j, j′] = val # push the element to the conjugate positions
            end
        end
        f = eigen(h)
        # save only target states
        energies = f.values[n_target_min:n_target_max]
        
        q = energies[n_w÷2+1] > (energies[n_w÷2] + energies[n_w÷2+2])/2 # true if the edge state branch is above the mean value
        
        # Lower band
        c = f.vectors[:, n_target_min:(n_target_min + n_w÷2 + !q - 1)]
        n_levels = size(c, 2)
        x = Matrix{Float64}(undef, n_levels, n_levels)
        for n in 1:n_levels
            for n′ in n:n_levels
                x[n′, n] = x[n, n′] = N*sum(c[j, n] * (π/2 * c[j, n′] - 8/π * sum(c[j′, n′]*j*j′/(j^2-j′^2)^2 for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
            end
        end
        pos_lower[z], d = eigen(x)
        ε_lower[z] = [dˣ.^2 ⋅ energies[1:(n_w÷2 + !q)] for dˣ in eachcol(d)]
        for i in 1:n_levels
            wf_lower[z][:, i] = abs2.(sum(d[j, i] * make_sine_state(coords, c[:, j]; n=N) for j = 1:n_levels))
        end

        # Higher band
        c = f.vectors[:, (n_target_min + n_w÷2 + !q):n_target_max]
        n_levels = size(c, 2)
        x = Matrix{Float64}(undef, n_levels, n_levels)
        for n in 1:n_levels
            for n′ in n:n_levels
                x[n′, n] = x[n, n′] = N*sum(c[j, n] * (π/2 * c[j, n′] - 8/π * sum(c[j′, n′]*j*j′/(j^2-j′^2)^2 for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
            end
        end
        pos_higher[z], d = eigen(x)
        ε_higher[z] = [dˣ.^2 ⋅ energies[(n_w÷2 + !q+1):end] for dˣ in eachcol(d)]
        for i in 1:n_levels
            wf_higher[z][:, i] = abs2.(sum(d[j, i] * make_sine_state(coords, c[:, j]; n=N) for j = 1:n_levels))
        end
    end
    return pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher
end

"Reconstruct the periodic coordinate space wavefunction 𝜓(𝑥) = ∑ⱼ𝑐ⱼexp(i2𝑗𝑥/𝑛) / √𝑛π"
function make_exp_state(x::AbstractVector{<:Real}, coeffs::AbstractVector{<:Number}; n)
    ψ = zeros(eltype(coeffs), length(x))
    n_j = (length(coeffs) - 1) ÷ 2
    for j in -n_j:n_j
        @. ψ += coeffs[j+n_j+1] * cis(2j/n * x)
    end
    return ψ ./ sqrt(n*π)
end

"Reconstruct the non-periodic coordinate space wavefunction 𝜓(𝑥) = ∑ⱼ𝑐ⱼsin(𝑗𝑥/𝑛) / √(𝑛π/2)"
function make_sine_state(x::AbstractVector{<:Real}, coeffs::AbstractVector{<:Number}; n)
    ψ = zeros(eltype(coeffs), length(x))
    for (j, c) in enumerate(coeffs)
        @. ψ += c * sin(j/n * x)
    end
    return ψ ./ sqrt(n*π/2)
end

"""

"""
function compute_wannier_centres_periodic(; N::Integer, n_max::Integer, n_target::Integer, phases::AbstractVector{<:Real}, gₗ::Real, Vₗ::Real)
    n_w = 2N # number of Wannier levels
    n_target_min = (n_target-1) * 2N + 1
    n_target_max = n_target_min + 2N - 1

    n_j = n_max * 2N # number of indices 𝑗 to use for constructing the unperturbed Hamiltonian

    # h = BM.BandedMatrix(BM.Zeros{ComplexF64}(2n_j + 1, 2n_j + 1), (2N, 2N))   ### Here and below we comment out the sparse version which sometimes gives LAPACKException(22)
    # h[BM.band(0)] .= [(2j/N)^2 + (gₗ + Vₗ)/2 for j = -n_j:n_j]
    # h[BM.band(-2N)] .= h[BM.band(2N)] .= gₗ/4
    h = zeros(ComplexF64, 2n_j + 1, 2n_j + 1)
    h = diagm(0 => ComplexF64[(2j/N)^2 + (gₗ + Vₗ)/2 for j = -n_j:n_j])
    h[diagind(h, -2N)] .= h[diagind(h, 2N)] .= gₗ/4

    energies = Matrix{Float64}(undef, n_w, length(phases))
    c = Matrix{ComplexF64}(undef, 2n_j+1, N) # eigenvectors of ℎ

    coords = range(0, N*pi, length=50N) # x's for wavefunctions
    pos_lower = Matrix{Float64}(undef, N, length(phases))
    pos_higher = Matrix{Float64}(undef, N, length(phases))
    ε_lower = Matrix{Float64}(undef, N, length(phases))
    ε_higher = Matrix{Float64}(undef, N, length(phases))
    wf_lower = Array{Float64, 3}(undef, length(coords), N, length(phases))
    wf_higher = Array{Float64, 3}(undef, length(coords), N, length(phases))

    x = Matrix{ComplexF64}(undef, N, N) # position operator
    d = Matrix{ComplexF64}(undef, N, N) # matrix of eigenvectors of the position operator
    pos_complex = Vector{Float64}(undef, N) # eigenvalues of the position operator; we will be taking their angles
    
    for (z, ϕ) in enumerate(phases)
        # h[BM.band(-N)] .= Vₗ/4 * cis(2ϕ)
        # h[BM.band(N)]  .= Vₗ/4 * cis(-2ϕ)
        h[diagind(h, -N)] .= Vₗ/4 * cis(2ϕ)
        h[diagind(h, N)]  .= Vₗ/4 * cis(-2ϕ)
        f = eigen(h)
        energies[:, z] = f.values[n_target_min:n_target_max]
        
        # Lower band
        # c = view(f[2], n_target_min:n_target_min + N - 1)
        c .= f.vectors[:, n_target_min:n_target_min + N - 1]
        for n in 1:N
            for n′ in 1:N
                # x[n′, n] = sum(c[n′][j+1]' * c[n][j] for j = 1:2n_j)
                x[n′, n] = sum(c[j+1, n′]' * c[j, n] for j = 1:2n_j)
            end
        end
        _, d, pos_complex = schur(x)
        pos_lower[:, z] = sort(@. (angle(pos_complex) + π) / 2π * N*π)
        ε_lower[:, z] = [abs2.(dˣ) ⋅ energies[1:N, z] for dˣ in eachcol(d)]
        for i in 1:N
            wf_lower[:, i, z] = abs2.(sum(d[j, i] * make_exp_state(coords, c[:, j]; n=N) for j = 1:N))
        end

        # Higher band
        # c = view(f[2], :, n_target_min + N:n_target_max)
        c .= f.vectors[:, n_target_min + N:n_target_max]
        for n in 1:N
            for n′ in 1:N
                # x[n′, n] = sum(c[n′][j+1]' * c[n][j] for j = 1:2n_j)
                x[n′, n] = sum(c[j+1, n′]' * c[j, n] for j = 1:2n_j)
            end
        end
        _, d, pos_complex = schur(x)
        pos_higher[:, z] = sort(@. (angle(pos_complex) + π) / 2π * N*π)
        ε_higher[:, z] = [abs2.(dˣ) ⋅ energies[N+1:2N, z] for dˣ in eachcol(d)]
        for i in 1:N
            wf_higher[:, i, z] = abs2.(sum(d[j, i] * make_exp_state(coords, c[:, j]; n=N) for j = 1:N))
        end
    end
    return energies, pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher
end

"""

"""
function compute_wannier_centres_qc(; n_levels::Integer, phases::AbstractVector{<:Real}, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real, χₗ::Real, χₛ::Real)
    X(j′, j, s) = 16s*j*j′ / (π*((j-j′)^2-4s^2)*((j+j′)^2-4s^2))
    
    n_j = 5n_levels # number of indices 𝑗 to use for constructing the Hamiltonian
    H = zeros(n_j, n_j)

    n_w = 2s - 1 # number of Wannier levels
    pos_lower = [Float64[] for _ in 1:length(phases)]
    pos_higher = [Float64[] for _ in 1:length(phases)]
    ε_lower = [Float64[] for _ in 1:length(phases)]
    ε_higher = [Float64[] for _ in 1:length(phases)]

    x = zeros(n_w÷2 + 1, n_w÷2 + 1)
    
    for (z, ϕ) in enumerate(phases)
        for j in 1:n_j
            for j′ in j:n_j
                val = 0.0
                if isodd(j′ + j)
                    val += λₗAₗ*X(j′, j, s)*sin(χₗ - ϕ) + λₛAₛ*X(j′, j, 2s)*sin(χₛ)
                else
                    # check diagonals "\"
                    if j′ == j
                        val += j^2 / 8M
                    elseif j′ == j - 2s || j′ == j + 2s
                        val += λₗAₗ * cos(χₗ - ϕ) / 2
                    elseif j′ == j - 4s || j′ == j + 4s
                        val += λₛAₛ * cos(χₛ) / 2
                    end
                    # check anti-diagonals "/"
                    if j′ == -j + 2s
                        val += -λₗAₗ * cos(χₗ - ϕ) / 2
                    elseif j′ == -j + 4s
                        val += -λₛAₛ * cos(χₛ) / 2
                    end
                end
                H[j′, j] = H[j, j′] = val # push the element to the conjugate positions
            end
        end

        f = eigen(H, sortby=-) # sort in descending order
        # save only target states
        energies = f.values[1:n_w]
        
        q = energies[n_w÷2+1] < (energies[n_w÷2] + energies[n_w÷2+2])/2 # true if the edge state branch is below the mean value

        # Lower band
        c = f.vectors[:, 1:(n_w÷2 + !q)]
        n_levels = size(c, 2)
        x = Matrix{Float64}(undef, n_levels, n_levels)
        for n in 1:n_levels
            for n′ in n:n_levels
                x[n′, n] = x[n, n′] = 2*sum(c[j, n] * (π/2 * c[j, n′] - 8/π * sum(c[j′, n′]*j*j′/(j^2-j′^2)^2 for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
            end
        end
        pos_lower[z], d = eigen(x)
        ε_lower[z] = [dˣ.^2 ⋅ energies[1:(n_w÷2 + !q)] for dˣ in eachcol(d)]

        # Higher band
        c = f.vectors[:, (n_w÷2 + !q + 1):n_w]
        n_levels = size(c, 2)
        x = Matrix{Float64}(undef, n_levels, n_levels)
        for n in 1:n_levels
            for n′ in n:n_levels
                x[n′, n] = x[n, n′] = 2*sum(c[j, n] * (π/2 * c[j, n′] - 8/π * sum(c[j′, n′]*j*j′/(j^2-j′^2)^2 for j′ = (iseven(j) ? 1 : 2):2:n_j)) for j = 1:n_j)
            end
        end
        pos_higher[z], d = eigen(x)
        ε_higher[z] = [dˣ.^2 ⋅ energies[(n_w÷2 + !q + 1):end] for dˣ in eachcol(d)]
    end
    return pos_lower, pos_higher, ε_lower, ε_higher
end

"""

"""
function compute_wannier_centres_qc_periodic(; phases::AbstractVector{<:Real}, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real, χₛ::Real, χₗ::Real)
    n_j = 10s

    h = BM.BandedMatrix(BM.Zeros{ComplexF64}(2n_j + 1, 2n_j + 1), (2s, 2s))
    h[BM.band(0)] .= [j^2 / M for j = -n_j:n_j]
    h[BM.band(-2s)] .= λₛAₛ * cis(χₛ)
    h[BM.band(2s)]  .= λₛAₛ * cis(-χₛ)

    energies = Matrix{Float64}(undef, 2s, length(phases))

    θ = range(0, 2π, length=40s) # x's for wavefunctions
    pos_lower = Matrix{Float64}(undef, s, length(phases))
    pos_higher = Matrix{Float64}(undef, s, length(phases))
    ε_lower = Matrix{Float64}(undef, s, length(phases))
    ε_higher = Matrix{Float64}(undef, s, length(phases))
    wf_lower = Array{Float64,3}(undef, length(θ), s, length(phases))
    wf_higher = Array{Float64,3}(undef, length(θ), s, length(phases))

    x = Matrix{ComplexF64}(undef, s, s) # position operator
    d = Matrix{ComplexF64}(undef, s, s) # matrix of eigenvectors of the position operator
    pos_complex = Vector{Float64}(undef, s) # eigenvalues of the position operator; we will be taking their angles
    
    for (z, ϕ) in enumerate(phases)
        h[BM.band(-s)] .= λₗAₗ * cis( χₗ - ϕ)
        h[BM.band(s)]  .= λₗAₗ * cis(-χₗ + ϕ)
        f = eigsolve(h, 2s, :LR; krylovdim=2n_j+1)
        energies[:, z] = f[1][1:2s]
        
        # Higher band
        c = view(f[2], 1:2s)
        for n in 1:s
            for n′ in 1:s
                x[n′, n] = sum(c[n′][j+1]' * c[n][j] for j = 1:2n_j) 
            end
        end
        _, d, pos_complex = schur(x)
        pos_real = angle.(pos_complex) .+ π
        sp = sortperm(pos_real)
        pos_higher[:, z] = pos_real[sp]
        Base.permutecols!!(d, sp)       # sort the eigenvalues in the same way
        ε_higher[:, z] = [abs2.(dˣ) ⋅ energies[1:s, z] for dˣ in eachcol(d)]
        for i in 1:s
            wf_higher[:, i, z] = abs2.(sum(d[j, i] * make_exp_state(θ, c[j]; n=2) for j = 1:s))
        end

        # Lower band
        c = view(f[2], 1+s:2s)
        for n in 1:s
            for n′ in 1:s
                x[n′, n] = sum(c[n′][j+1]' * c[n][j] for j = 1:2n_j) 
            end
        end
        _, d, pos_complex = schur(x)
        pos_real = angle.(pos_complex) .+ π
        sp = sortperm(pos_real)
        pos_lower[:, z] = pos_real[sp]
        Base.permutecols!!(d, sp)       # sort the eigenvalues in the same way
        ε_lower[:, z] = [abs2.(dˣ) ⋅ energies[s+1:2s, z] for dˣ in eachcol(d)]
        for i in 1:s
            wf_lower[:, i, z] = abs2.(sum(d[j, i] * make_exp_state(θ, c[j]; n=2) for j = 1:s))
        end
    end
    return energies, pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher
end