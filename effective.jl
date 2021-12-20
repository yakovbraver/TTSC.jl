using UnPack
using BandedMatrices: BandedMatrix, band
using KrylovKit: eigsolve

"""
Hamiltonian (S17) defined here as
    𝐻 = (ℎ₀ - Ω𝐼₀) + (𝐼 - 𝐼₀)²/2𝑀 + λL⋅cos(𝑠ϑ + χL) + λS⋅cos(2𝑠ϑ + χS)
"""
Base.@kwdef mutable struct EffectiveHamiltonian
    h0::Float64
    Ω::Float64
    I0::Float64
    M::Float64
    λL::Float64
    χL::Float64
    λS::Float64
    χS::Float64
    s::Float64
end

"Calculate energies for the values of the action-angle variables given in vectors `I` and `ϑ`"
function energies_in_phasespace(H::EffectiveHamiltonian, ϑ::AbstractVector, I::AbstractVector)
    @unpack h0, Ω, I0, M, λL, χL, λS, χS, s = H
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = h0 - Ω*I0 + (I[i]-I0)^2/2M + λL*cos(s*ϑ[t] + χL) + λS*cos(2s*ϑ[t] + χS)
    end
    E
end

"""
Calculate `n_bands` energy bands of Hamiltonian (S32) sweeping over the adiabatic `phases`.
In the returned Matrix of bands, columns numerate the adiabatic phase, while rows numerate eigenvalues.
Rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
Rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case 𝑘 = 1.
"""
function get_bands(; n_bands::Integer, phases::AbstractVector, M::Real, λL::Real, λS::Real)
    n_j = 4n_bands # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
    # Hamiltonian matrix
    H = BandedMatrix{ComplexF64}(undef, (2n_j + 1, 2n_j + 1), (2, 2))
    H[band(-2)] .= fill(λS, 2n_j-1)
    H[band(2)]  .= fill(λS, 2n_j-1)
    
    bands = Matrix{Float64}(undef, 2n_bands, length(phases))

    for k in [0, 1]
        H[band(0)] .= [(2j + k)^2 / M for j = -n_j:n_j]
        # `a` and `b` control where to place the eigenvalues depedning on `k`; see description of `bands`
        a = k*n_bands + 1 
        b = a+n_bands - 1
        for (i, ϕ) in enumerate(phases)
            H[band(-1)] .= fill(λL*cis(-ϕ), 2n_j)
            H[band(1)]  .= fill(λL*cis(ϕ), 2n_j)
            vals, _, _ = eigsolve(H, n_bands, :SR)
            bands[a:b, i] .= vals[1:n_bands]
        end
    end
    bands / 4π # we do not include the factor 4π in the diagonalisation problem and restore it here
end