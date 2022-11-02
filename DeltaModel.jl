module DeltaModel

import Roots
using LinearAlgebra: eigen, schur, ⋅, diagm, diagind, eigvals

"""
A type representing the unperturbed Hamiltonian
    ℎ = 𝑝² + 𝜆 ∑ₙ𝛿(𝑥 - 𝑛𝑎/3) + 𝑈 ∑ₙ𝑔ₙ(𝑥)cos(𝜑ₓ + 2π𝑛/3).
"""
mutable struct UnperturbedHamiltonian
    N::Int # number of lattice cells
    a::Float64
    λ::Float64
    U::Float64
    isperiodic::Bool
    φₓ::Vector{Float64}
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue of the band of interest at `j`th phase, `i` ∈ [1, `N`], `j` ∈ [1, `length(φₓ)`]
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector of the band of interest at `j`th phase
end

"""
Construct an `UnperturbedHamiltonian` object.
"""
function UnperturbedHamiltonian(n_cells::Integer; a::Real, λ::Real, U::Real, isperiodic::Bool, φₓ::AbstractVector{<:Real})
    E = Matrix{Float64}(undef, n_cells, length(φₓ))
    c = Array{ComplexF64,3}(undef, 6, n_cells, length(φₓ))

    UnperturbedHamiltonian(Int(n_cells), Float64(a), Float64(λ), Float64(U), isperiodic, collect(Float64, φₓ), E, c)
end

"Return 𝑔ₙ(𝑥) which realises the pumping protocol."
function 𝑔(x; n, a)
    Int( n/3 <= (x % a)/a < (n+1)/3 )
end

"Return cos(𝑘𝑎)/2𝜅₁𝜅₂𝜅₃ for a given energy `ε` and phase `φ`."
function cos_ka(ε::Real; φ::Real, uh::UnperturbedHamiltonian)
    (;a, U, λ) = uh
    κ = [sqrt(ε - U*cos(φ + 2π*n/3)) for n = 0:2]
    s = sin.(κ*a/3)
    c = cos.(κ*a/3)
    return (-κ[1]^2*s[1] * (s[2]*(λ*s[3]+κ[3]c[3]) + κ[2]s[3]c[2]) +
            κ[1]c[1] * (s[2] * (3λ*κ[3]c[3] - s[3]*(κ[2]^2+κ[3]^2-2λ^2)) + κ[2]c[2]*(3λ*s[3]+2κ[3]c[3])) +
            s[1] * (s[2]*(λ*s[3]*(-κ[2]^2-κ[3]^2+λ^2) + κ[3]c[3]*(2λ^2-κ[2]^2)) + κ[2]c[2]*(3λ*κ[3]c[3] - s[3]*(κ[3]^2-2λ^2)))
    ) / 2κ[1]κ[2]κ[3]
end

"Find allowed energies for the Hamiltonian `uh` at each phase for a band bracketed in energy by `bounds`."
function diagonalise!(uh::UnperturbedHamiltonian, bounds::Tuple{<:Real, <:Real})
    (;N) = uh
    for (i, φ) in enumerate(uh.φₓ)
        mid = N÷2 + 1
        for n in 0:mid
            uh.E[n+1, i] = Roots.find_zero(ε -> cos_ka(ε; φ, uh) - cos(2π*n/N), bounds, Roots.A42(), rtol=1e-5)
        end
        uh.E[mid+1:N, i] = uh.E[N-mid+1:-1:2, i]
    end
end

end