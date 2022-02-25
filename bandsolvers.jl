import BandedMatrices as BM
using SparseArrays: sparse
using KrylovKit: eigsolve

"""
Calculate `n_bands` energy bands of Hamiltonian (S32) sweeping over the adiabatic `phases` (φₜ in (S32)).
In the returned matrix of bands, columns enumerate the adiabatic phases, while rows enumerate eigenvalues.
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

function compute_bands_with_boundary(; n_bands::Integer, phases::AbstractVector, M::Real, λₗAₗ::Real, λₛAₛ::Real)
    n_j = 20 # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
    X(j′, j) = 32*j*j′ / (π*((j-j′)^2-16)*((j+j′)^2-16))

    # # Hamiltonian matrix
    # H_rows, H_cols, H_vals = Float64[], Float64[], Float64[]
    # sizehint!(H_rows, (2n_j-1)^2)    
    # sizehint!(H_cols, (2n_j-1)^2)    
    # sizehint!(H_vals, (2n_j-1)^2)
    # # Fill X's; skip zeroth row and column
    # for j′ = [-n_j:1, 1:n_j]
    #     for j = [-n_j+!(j′&1):-1, 1+!(j′&1):n_j] # start with `-n_j+1` if `j′` is even, and just `-n_j` otherwise
    #         push!(H_cols, j′); push!(H_rows, j)
    #         push!(H_vals, X(j′, j)) * sin(ϕ)
    #     end
    # end
    # # Fill the main diagonal, skip cells (±4, ±4) and (±2, ±2)
    # for j = [-n_j:-5, 3, -1:1, 3, 5:n_j]
    #     push!(H_cols, j); push!(H_rows, j)
    #     push!(H_vals, j^2/8M)
    # end
    # # Fill (±4, ±4)
    # push!(H_cols, 4); push!(H_rows, 4)
    # push!(H_vals, 4^2/8M - λₛAₛ/2)
    # push!(H_cols, -4); push!(H_rows, -4)
    # push!(H_vals, 4^2/8M - λₛAₛ/2)
    # # Fill (±2, ±2)
    # push!(H_cols, 2); push!(H_rows, 2)
    # push!(H_vals, 2^2/8M - λₗAₗ*cos(ϕ)/2)
    # push!(H_cols, -2); push!(H_rows, -2)
    # push!(H_vals, 2^2/8M - λₗAₗ*cos(ϕ)/2)
    # # Fill the main diagonal, skip cells (±4, ±4) and (±2, ±2)
    # for j = [-n_j:-5, 3, -1:1, 3, 5:n_j]
    #     push!(H_cols, j); push!(H_rows, j)
    #     push!(H_vals, j^2/8M)
    # end
    # # Fill λₛAₛ diagonals
    # for j′ = [-n_j:-8, -7, -5:-3, -2:-1, 1:n_j-8]
    #     j = j′ + 8
    #     push!(H_cols, j′); push!(H_rows, j)
    #     push!(H_vals, λₛAₛ/2)
    #     # place conjugate element
    #     push!(H_cols, j); push!(H_rows, j′)
    #     push!(H_vals, λₛAₛ/2)
    # end
    # # Fill (2, -6) and (-6, 2)
    # push!(H_rows, 2); push!(H_cols, -6); 
    # push!(H_vals, λₛAₛ/2 - λₗAₗ*cos(ϕ)/2)
    # push!(H_cols, -6); push!(H_rows, 2)
    # push!(H_vals, λₛAₛ/2 - λₗAₗ*cos(ϕ)/2)
    # # Fill (-2, 6) and (6, -2)
    # push!(H_rows, -2); push!(H_cols, 6); 
    # push!(H_vals, λₛAₛ/2 - λₗAₗ*cos(ϕ)/2)
    # push!(H_cols, 6); push!(H_rows, -2)
    # push!(H_vals, λₛAₛ/2 - λₗAₗ*cos(ϕ)/2)

    H_rows, H_cols, H_vals = Float64[], Float64[], Float64[]
    sizehint!(H_rows, (2n_j-1)^2)    
    sizehint!(H_cols, (2n_j-1)^2)    
    sizehint!(H_vals, (2n_j-1)^2)
    for j in [-n_j:-1; 1:n_j]
        cols = j < 0 ? [j:-1; 1:n_j] : j:n_j
        for j′ in cols
            val = 0.0
            if abs(j′ + j) % 2 == 1 # if `j′ + j` is odd
                val += X(j′, j) * sin(0)
            else
                # check diagonals "\"
                if j′ == j
                    val += j^2 / 8M
                elseif j′ == j - 4 || j′ == j + 4
                    val += λₗAₗ * cos(0) / 2
                elseif j′ == j - 8 || j′ == j + 8
                    val += λₛAₛ / 2
                end
                # check diagonals "/"
                if j′ == -j - 4 || j′ == -j + 4
                    val += -λₗAₗ * cos(0) / 2
                elseif j′ == -j - 8 || j′ == -j + 8
                    val += -λₛAₛ / 2
                end
            end
            if val != 0 || abs(j′ + j) % 2 == 1
                # Push the non-zero element to the conjugate positions; also "shift" the indices because they must start from 1
                push!(H_rows, j′+n_j+1); push!(H_cols, j+n_j+1)
                push!(H_vals, val)
                # println("($j′, $j): $val")
                if j′ != j
                    push!(H_rows, j+n_j+1); push!(H_cols, j′+n_j+1)
                    push!(H_vals, val)
                end
            end
        end
    end
    # return sparse(H_rows, H_cols, H_vals)

    H = sparse(H_rows, H_cols, H_vals)
    bands = Matrix{Float64}(undef, n_bands, length(phases))
    # format of eigenstates is `states[phase][eigenstate_number][j]`
    states = [[Vector{Float64}(undef, 2n_j+1) for _ in 1:n_bands] for _ in 1:length(phases)]  # eigenvectors of 𝐻, we will save `n_bands` of them, and each will have `2n_j+1` components
    for (i, ϕ) in enumerate(phases)
        vals_counter = 1
        for j in [-n_j:-1; 1:n_j]
            cols = j < 0 ? [j:-1; 1:n_j] : j:n_j
            for j′ in cols
                val = 0.0
                if abs(j′ + j) % 2 == 1 # if `j′ + j` is odd
                    val += X(j′, j) * sin(ϕ)
                else
                    # check diagonals "\"
                    if j′ == j
                        val += j^2 / 8M
                    elseif j′ == j - 4 || j′ == j + 4
                        val += λₗAₗ * cos(ϕ) / 2
                    elseif j′ == j - 8 || j′ == j + 8
                        val += λₛAₛ / 2
                    end
                    # check diagonals "/"
                    if j′ == -j - 4 || j′ == -j + 4
                        val += -λₗAₗ * cos(ϕ) / 2
                    elseif j′ == -j - 8 || j′ == -j + 8
                        val += -λₛAₛ / 2
                    end
                end
                if val != 0 || abs(j′ + j) % 2 == 1
                    H_vals[vals_counter] = val
                    vals_counter += 1
                    if j′ != j
                        H_vals[vals_counter] = val
                        vals_counter += 1
                    end
                end
            end
        end
        H = sparse(H_rows, H_cols, H_vals)
        vals, vecs, info = eigsolve(H, n_bands, :LR; krylovdim=30)
        if info.converged < n_bands
            @warn "Only $(info.converged) eigenvalues out of $(n_bands) converged when diagonalising 𝐻. "*
                  "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
        end
        states[i] .= vecs[1:n_bands]
        bands[:, i] .= vals[1:n_bands]
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
function compute_floquet_bands(; n_min::Integer, n_max::Integer, phases::AbstractVector, s::Integer, l::Real, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real, pumptype::Symbol)
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
    
    Eₖ = Matrix{Float64}(undef, 2Δn, length(phases)) # eigenvalues of 𝐻ₖ (Floquet quasi-energies) that will be saved; size is twice `Δn` for the two values of 𝑘
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
        a_Hₖ = (k > 0)*Δn + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_Hₖ = a_Hₖ+Δn - 1
        for (z, ϕ) in enumerate(phases)
            if pumptype != :time || z == 1 # If pupming is not time-only, ℎₖ has to be diagonalised on each iteration. If it's time-only, then we diagonalise only once, at `z == 1`.
                hₖ[BM.band(-1)] .= Vₗ/4 * cis(2ϕ)
                hₖ[BM.band(1)]  .= Vₗ/4 * cis(-2ϕ)
                vals, vecs, info = eigsolve(hₖ, 2n_max, :SR; krylovdim=2n_j+1)
                if info.converged < 2n_max
                    @warn "Only $(info.converged) eigenvalues out of $(2n_max) converged when diagonalising ℎₖ. "*
                          "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
                end
                # save only energies and states for levels from `2n_min` to `2n_max`
                ϵₖ[a_hₖ:b_hₖ, z] = vals[2n_min-1:2n_max]
                cₖ .= vecs[2n_min-1:2n_max]
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
            vals, _, info = eigsolve(Hₖ, Δn, :LR; krylovdim=Hₖ_dim)
            if info.converged < Δn
                @warn "Only $(info.converged) eigenvalues out of $(Δn) converged when diagonalising 𝐻ₖ. "*
                      "Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
            end
            Eₖ[a_Hₖ:b_Hₖ, z] .= vals[1:Δn]
        end
    end
    return ϵₖ, Eₖ
end