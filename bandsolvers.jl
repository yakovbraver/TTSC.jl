import BandedMatrices as BM
using SparseArrays: sparse
using KrylovKit: eigsolve

"""
Calculate `n_bands` energy bands of Hamiltonian (S32) sweeping over the adiabatic `phases` (φₜ in (S32)).
In the returned matrix of bands, columns numerate the adiabatic phase, while rows numerate eigenvalues.
Rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
Rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case λₗAₗcos(sϑ+φₜ) leads to 𝑘 = s/2.
"""
function compute_secular_bands(; n_bands::Integer, phases::AbstractVector, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real)
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
            H[BM.band(-1)] .= λₗAₗ*cis(-ϕ)
            H[BM.band(1)]  .= λₗAₗ*cis(ϕ)
            vals, _, _ = eigsolve(H, n_bands, :LR; krylovdim=n_bands+10)
            bands[a:b, i] .= vals[1:n_bands]
        end
    end
    return bands / 2 # restore the omitted factor
end

"""
Calculate `n_bands` energy bands of the Floquet Hamiltonian (S20) sweeping over the adiabatic `phases` φₓ and φₜ.
Return a tuple of a matrix `ϵₖ` of `2n_bands` bands of ℎₖ and a matrix `Eₖ` of `n_bands` bands of 𝐻ₖ.
In the returned matrices, columns numerate the adiabatic phases, while rows numerate eigenvalues.
In `Eₖ`, rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
In `Eₖ`, rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case Vₗcos²(x+φₓ) leads to 𝑘 = 2/2 = 1.
The dimension of the constructed 𝐻ₖ matrix will be `2n_bands`, hence that many eigenvalues of ℎₖ will be required. This in turn
requires constructing ℎₖ of dimension `4n_bands`.
Type of pumping is controlled via `pumptype`: `:time` for temporal, `:space` for spatial, or anything else for simultaneous space-time pumping.
Note that if `pumptype==:time`, ℎₖ is diagonalised only once (as the spatial phase is constant), hence only the first column of `ϵₖ` are populated.
"""
function compute_floquet_bands(; n_bands::Integer, phases::AbstractVector, s::Integer, l::Real, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
    n_j = 2n_bands  # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
    hₖ = BM.BandedMatrix(BM.Zeros{ComplexF64}(2n_j+1, 2n_j+1), (2l, 2l))
    # fill the off-diagonals with binomial numbers; the diagonal is treated in the `k` loop
    for n in 1:l
        hₖ[BM.band(2n)] .= hₖ[BM.band(-2n)] .= gₗ / 4^l * binomial(2l, l-n)
    end
    
    ϵₖ = Matrix{Float64}(undef, 2n_j, length(phases)) # eigenvalues of ℎₖ (eigenenergies of the unperturbed Hamiltonian)
    cₖ = [Vector{ComplexF64}(undef, 2n_j+1) for _ in 1:2n_bands]  # eigenvectors of ℎₖ, we will save `2n_bands` of them, and each will have `2n_j+1` components
    
    Eₖ = Matrix{Float64}(undef, 2n_bands, length(phases)) # eigenvalues of 𝐻ₖ (Floquet quasi-energies) that will be saved; size is twice `n_bands` for the two values of `k`
    Hₖ_dim = 2n_bands # dimension of the constructed 𝐻ₖ matrix (twice larger than the number of requested quasi-energies)
    n_Hₖ_nonzeros = 9Hₖ_dim - 24s # number of non-zero elements in 𝐻ₖ
    Hₖ_rows = Vector{Int}(undef, n_Hₖ_nonzeros)
    Hₖ_cols = Vector{Int}(undef, n_Hₖ_nonzeros)
    Hₖ_vals = Vector{ComplexF64}(undef, n_Hₖ_nonzeros)
    for k in [0, 1] # iterate over the centre of BZ and then the boundary
        hₖ[BM.band(0)] .= [(2j + k)^2 + Vₗ/2 + gₗ / 4^l * binomial(2l, l) for j = -n_j:n_j]
        # `a_*` and `b_*` control where to place the eigenvalues of 𝐻ₖ and ℎₖ depedning on `k`; see function docstring
        a_Hₖ = (k > 0)*n_bands + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_Hₖ = a_Hₖ+n_bands - 1
        a_hₖ = (k > 0)*n_j + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_hₖ = a_hₖ+n_j - 1
        for (z, ϕ) in enumerate(phases)
            if pumptype != :time || z == 1 # If pupming is not time-only, ℎₖ has to be diagonalised on each iteration. If it's time-only, then we diagonalise only once, at `z == 1`.
                hₖ[BM.band(-1)] .= Vₗ/4 * cis(2ϕ)
                hₖ[BM.band(1)]  .= Vₗ/4 * cis(-2ϕ)
                vals, vecs, info = eigsolve(hₖ, n_j, :SR; krylovdim=n_j+10)
                if info.converged < n_j
                    @warn "Only $(info.converged) eigenvalues out of $(n_j) converged when diagonalising ℎₖ."*
                          "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
                end
                ϵₖ[a_hₖ:b_hₖ, z] = vals[1:n_j]
                cₖ .= vecs[1:n_j]
            end

            # Construct 𝐻ₖ
            p = 1 # a counter for placing elements to the vectors `Hₖ_*`
            for m in 1:Hₖ_dim
                # place the diagonal element (S25)
                Hₖ_rows[p] = Hₖ_cols[p] = m
                q = (pumptype == :time ? 1 : z) # If pumping is time-only, `ϵₖ[m, z]` is only calculated for `z == 1` (during diagonalisation of ℎₖ)
                Hₖ_vals[p] = ϵₖ[m, q] - ceil(m/2)*ω/s
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
                        Hₖ_vals[p] = (pumptype == :space ? λₗ/2 * j_sum : λₗ/2 * j_sum * cis(-ϕ)) # a check for space or space-time pumping
                    elseif pumptype == :time 
                        Hₖ_vals[p] *= cis(-(phases[2]-phases[1]))
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
            vals, _, info = eigsolve(Hₖ, n_bands, :LR; krylovdim=n_bands+10)
            if info.converged < n_bands
                @warn "Only $(info.converged) eigenvalues out of $(n_bands) converged when diagonalising 𝐻ₖ."*
                      "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
            end
            Eₖ[a_Hₖ:b_Hₖ, z] .= vals[1:n_bands]
        end
    end
    return ϵₖ, Eₖ
end