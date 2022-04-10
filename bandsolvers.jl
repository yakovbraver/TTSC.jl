import BandedMatrices as BM
using SparseArrays: sparse
using KrylovKit: eigsolve
using LinearAlgebra: eigen

"""
Calculate `n_bands` of "quasiclassical" energy bands of Hamiltonian (S32) sweeping over the adiabatic `phases` (φₜ in (S32)).
In the returned matrix of bands, columns enumerate the adiabatic phases, while rows enumerate eigenvalues.
Rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
Rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case λₗAₗcos(sϑ+φₜ) leads to 𝑘 = s/2.
"""
function compute_qc_bands(; n_bands::Integer, phases::AbstractVector{<:Real}, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real)
    n_j = 2n_bands # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
    # Hamiltonian matrix
    H = BM.BandedMatrix{ComplexF64}(undef, (2n_j + 1, 2n_j + 1), (2, 2))
    H[BM.band(-2)] .= λₛAₛ
    H[BM.band(2)]  .= λₛAₛ
    
    bands = Matrix{Float64}(undef, 2n_bands, length(phases))
    for k in [0, s÷2] # iterate over the centre of BZ and then the boundary
        H[BM.band(0)] .= [(2j + k)^2 / M for j = -n_j:n_j]
        # `a` and `b` control where to place the eigenvalues depedning on `k`; see function docstring
        a = (k > 0)*n_bands + 1 
        b = a+n_bands - 1
        for (i, ϕ) in enumerate(phases)
            H[BM.band(-1)] .= λₗAₗ*cis(-ϕ-pi/2)
            H[BM.band(1)]  .= λₗAₗ*cis(ϕ+pi/2)
            vals, _, _ = eigsolve(H, n_bands, :LR; krylovdim=n_bands+10)
            bands[a:b, i] .= vals[1:n_bands]
        end
    end
    return bands / 2 # restore the omitted factor
end

"""
Calculate all "quasiclassical" energy bands of Hamiltonian (S32) with boundaries, sweeping over the adiabatic `phases` (φₜ in (S32)).
Return a tuple (`bands`, `states`): `bands[:, p]` stores eigenenergies at `p`th phase, while `states[p][:, m]` stores `m`th eigenvector at `p`th phase.
Bands and states are sorted in energy-descending order so that for `M` negative, the bands of interest will be the first ones.
Parameter `n` is the number of cells in the lattice; the eigenfunctions will be calculated as a superposition of functions sin(𝑗𝑥/𝑛) / √(𝑛π/2)
"""
function compute_qc_bands_with_boundary(; phases::AbstractVector{<:Real}, M::Real, λₗAₗ::Real, λₛAₛ::Real, n::Integer=2)    
    X(j′, j) = 16n*j*j′ / (π*((j-j′)^2-(2n)^2)*((j+j′)^2-(2n)^2))
    
    n_j = 300 # number of indices 𝑗 to use for constructing the Hamiltonian
    H = zeros(n_j, n_j)
    # for storing eigenstates and eigenvectors, see function docstring for format
    bands = Matrix{Float64}(undef, n_j, length(phases))
    states = [Matrix{Float64}(undef, n_j, n_j) for _ in 1:length(phases)]
    for (i, ϕ) in enumerate(phases)
        for j in 1:n_j
            for j′ in 1:n_j
                val = 0.0
                if abs(j′ + j) % 2 == 1 # if `j′ + j` is odd
                    val += λₗAₗ * X(j′, j) * sin(ϕ)
                else
                    # check diagonals "\"
                    if j′ == j
                        val += j^2 / (2M * n^2)
                    elseif j′ == j - 2n || j′ == j + 2n
                        val += λₗAₗ * cos(ϕ) / 2
                    elseif j′ == j - 4n || j′ == j + 4n
                        val += λₛAₛ / 2
                    end
                    # check diagonals "/"
                    if j′ == -j - 2n || j′ == -j + 2n
                        val += -λₗAₗ * cos(ϕ) / 2
                    elseif j′ == -j - 4n || j′ == -j + 4n
                        val += -λₛAₛ / 2
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
Note that if `pumptype==:time`, ℎₖ is diagonalised only once (as the spatial phase is constant), hence only the first column of `ϵₖ` is populated.
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
                    for p in 2:length(phases)
                        ϵₖ[a_hₖ:b_hₖ, p] = ϵₖ[a_hₖ:b_hₖ, 1]
                    end
                end
            end

            # Construct 𝐻ₖ
            p = 1 # a counter for placing elements to the vectors `Hₖ_*`
            for m in 1:Hₖ_dim
                # place the diagonal element (S25)
                Hₖ_rows[p] = Hₖ_cols[p] = m
                q = (pumptype == :time ? 1 : z) # If pumping is time-only, `ϵₖ[m, z]` is only calculated for `z == 1` (during diagonalisation of ℎₖ)
                Hₖ_vals[p] = ϵₖ[m+a_hₖ-1, q] - ceil((m+2n_min-2)/2)*ω/s
                p += 1

                # place the elements of the long lattice (S26)
                for i in 1:2
                    m′ = 2s + 2(ceil(Int, m/2)-1) + i
                    m′ > Hₖ_dim && break
                    Hₖ_rows[p] = m′
                    Hₖ_cols[p] = m
                    if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                        j_sum = sum( (cₖ[m′][j+2]/4 + cₖ[m′][j-2]/4 + cₖ[m′][j]/2)' * cₖ[m][j] for j = 3:2n_j-1 ) + 
                                     (cₖ[m′][3]/4 + cₖ[m′][1]/2)' * cₖ[m][1] +                # iteration j = 1
                                     (cₖ[m′][2n_j-1]/4 + cₖ[m′][2n_j+1]/2)' * cₖ[m][2n_j+1]   # iteration j = 2n_j+1
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
                        j_sum = sum( (-cₖ[m′][j+2]/4 - cₖ[m′][j-2]/4 + cₖ[m′][j]/2)' * cₖ[m][j] for j = 3:2n_j-1 ) + 
                                     (-cₖ[m′][3]/4 + cₖ[m′][1]/2)' * cₖ[m][1] +                # iteration j = 1
                                     (-cₖ[m′][2n_j-1]/4 + cₖ[m′][2n_j+1]/2)' * cₖ[m][2n_j+1]   # iteration j = 2n_j+1
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
are the Floquet quasienergies. Then, we sort then in descending order (as if we diagonalised the Hamiltonian) and find the permutation
that would undo this sorting. This permutation is applied to `E`.
The procedure yield fully correct results only if `E` has been calculated at zero perturbation. The perturbation may additionally change
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
WIP: Floquet spectrum with spatial boundaries.
Currently, calculation of the spatial spectrum is correct, but for the Floquet spectrum it is not.
"""
function compute_floquet_bands_with_boundary(;n::Integer, n_min::Integer, n_max::Integer, phases::AbstractVector{<:Real}, s::Integer, l::Real, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
    X(j′, j) = 16n*j*j′ / (π*((j-j′)^2-(2n)^2)*((j+j′)^2-(2n)^2))
    
    n_j = 600 # number of indices 𝑗 to use for constructing the Hamiltonian
    H = zeros(n_j, n_j)
    # G = 5 #n #+ 1
    G = repeat([5, 5, 5, 5, 5, 3, 3, 3], 7)
        push!(G, 3)
    Δn = n_max - n_min + 1
    
    # Eigenvalues of ℎₖ (eigenenergies of the unperturbed Hamiltonian).
    # We should store `2Δn` of them because each of the `Δn` levels are almost degenerate. To account for the two values of 𝑘, we use `4Δn`.
    ϵₖ = Matrix{Float64}(undef, Δn, length(phases))
    cₖ = Matrix{Float64}(undef, n_j, Δn)
    
    Eₖ = Matrix{Float64}(undef, Δn, length(phases)) # eigenvalues of 𝐻ₖ (Floquet quasi-energies) that will be saved; size is twice `Δn` for the two values of 𝑘
    Hₖ_dim = Δn # dimension of the constructed 𝐻ₖ matrix (twice larger than the number of requested quasi-energies)
   GG = 5
    n_Hₖ_nonzeros = 801# (4GG+1)*Hₖ_dim - 6GG^2*s # number of non-zero elements in 𝐻ₖ
   
    Hₖ_rows = zeros(Int, n_Hₖ_nonzeros)
    Hₖ_cols = zeros(Int, n_Hₖ_nonzeros)
    Hₖ_vals = zeros(ComplexF64, n_Hₖ_nonzeros)
    
    for (z, ϕ) in enumerate(phases)
        if pumptype != :time || z == 1 # If pupming is not time-only, ℎₖ has to be diagonalised on each iteration. If it's time-only, then we diagonalise only once, at `z == 1`.
            for j in 1:n_j
                for j′ in 1:n_j
                    val = 0.0
                    if abs(j′ + j) % 2 == 1 # if `j′ + j` is odd
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
                        # check diagonals "/"
                        if j′ == -j - 2n || j′ == -j + 2n
                            val += -Vₗ/2 * cos(2ϕ) / 2
                        elseif j′ == -j - 4n || j′ == -j + 4n
                            val += -gₗ/2 / 2
                        end
                    end
                    H[j′, j] = H[j, j′] = val # push the element to the conjugate positions
                end
            end
            f = eigen(H)
            # save only energies and states for levels from `2n_min` to `2n_max`
            ϵₖ[:, z] = f.values[n_min:n_max]
            cₖ .= f.vectors[:, n_min:n_max]
        end
        # return ϵₖ, cₖ
        # Construct 𝐻ₖ
        p = 1 # a counter for placing elements to the vectors `Hₖ_*`
        
        for m in 1:Hₖ_dim
            # place the diagonal element (S25)
            Hₖ_rows[p] = Hₖ_cols[p] = m
            q = (pumptype == :time ? 1 : z) # If pumping is time-only, `ϵₖ[m, z]` is only calculated for `z == 1` (during diagonalisation of ℎₖ)
            Hₖ_vals[p] = ϵₖ[m, q] - ceil(m/G[m])*ω/s
            p += 1

            # place the elements of the long lattice (S26)
            for i in 1:G[m]
                m′ = G[m]*s + G[m]*(ceil(Int, m/G[m])-1) + i
                m′ > Hₖ_dim && break
                Hₖ_rows[p] = m′
                Hₖ_cols[p] = m
                if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                    j_sum = sum( (cₖ[j+4n, m′]/4 + cₖ[j-4n, m′]/4 + cₖ[j, m′]/2) * cₖ[j, m] for j = 4n+1:n_j-4n ) + 
                            sum( (cₖ[j+4n, m′]/4 - cₖ[-j+4n, m′]/4 + cₖ[j, m′]/2) * cₖ[j, m] for j = 1:4n-1 ) +
                            (cₖ[4n+4n, m′]/4 + cₖ[4n, m′]/2) * cₖ[4n, m] + # iteration `j = 4n`
                            sum( (cₖ[j-4n, m′]/4 + cₖ[j, m′]/2) * cₖ[j, m] for j = n_j-4n+1:n_j )
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
            for i in 1:G[m]
                m′ = 2s*G[m] + G[m]*(ceil(Int, m/G[m])-1) + i
                m′ > Hₖ_dim && break
                Hₖ_rows[p] = m′
                Hₖ_cols[p] = m
                if pumptype != :time || z == 1 # If pumping is time-only, this may be calculated only once
                    j_sum = sum( (-cₖ[j+4n, m′]/4 - cₖ[j-4n, m′]/4 + cₖ[j, m′]/2) * cₖ[j, m] for j = 4n+1:n_j-4n ) + 
                            sum( (-cₖ[j+4n, m′]/4 + cₖ[-j+4n, m′]/4 + cₖ[j, m′]/2) * cₖ[j, m] for j = 1:4n-1) +
                            (-cₖ[4n+4n, m′]/4 + cₖ[4n, m′]/2) * cₖ[4n, m] + # iteration `j = 4n`
                            sum( (-cₖ[j-4n, m′]/4 + cₖ[j, m′]/2) * cₖ[j, m] for j = n_j-4n+1:n_j)
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
        # println(count(!=(0), Hₖ_rows))
        Hₖ = sparse(Hₖ_rows, Hₖ_cols, Hₖ_vals)
        vals, _, info = eigsolve(Hₖ, Δn, :LR; krylovdim=Hₖ_dim)
        if info.converged < Δn
            @warn "Only $(info.converged) eigenvalues out of $(Δn) converged when diagonalising 𝐻ₖ. "*
                  "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
        end
        Eₖ[:, z] .= vals[1:Δn]
    end
    return ϵₖ, Eₖ
end