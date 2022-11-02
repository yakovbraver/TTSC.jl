module DeltaModel

import Roots
using LinearAlgebra: eigen, schur, ⋅, diagm, diagind, eigvals, det

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

function system_matrix(ε::Real; φ::Real, ka::Real, uh::UnperturbedHamiltonian)
    (;a, U, λ) = uh
    κ = [sqrt(ε - U*cos(φ + 2π*n/3)) for n = 0:2]
    s = [sin(n*κ[m]*a/3) for n = 1:3, m = 1:3]
    c = [cos(n*κ[m]*a/3) for n = 1:3, m = 1:3]
    e = cis(-ka)
    [0                    -1                  0                    0                   s[3,3]e       c[3,3]e;
     s[1,1]               c[1,1]              -s[1,2]              -c[1,2]             0             0;
     0                    0                   s[2,2]               c[2,2]              -s[2,3]       -c[2,3];
     κ[1]                 -λ                  0                    0                   -κ[3]c[3,3]e  κ[3]s[3,3]e;
     -κ[1]c[1,1]-s[1,1]λ  κ[1]s[1,1]-c[1,1]λ  κ[2]c[1,2]           -κ[2]s[1,2]         0             0;
     0                    0                   -κ[2]c[2,2]-s[2,2]λ  κ[2]s[2,2]-c[2,2]λ  κ[3]c[2,3]    -κ[3]s[2,3]]
end

"Find allowed energies for the Hamiltonian `uh` at each phase for a band bracketed in energy by `bounds`."
function diagonalise!(uh::UnperturbedHamiltonian, bounds::Tuple{<:Real, <:Real})
    (;N, a, U, φₓ, E) = uh
    for (j, φ) in enumerate(φₓ)
        for i = 1:N
            ka = 2π*(i-1)/N
            # eigenenergies
            if i <= N÷2 + 1
                E[i, j] = Roots.find_zero(ε -> cos_ka(ε; φ, uh) - cos(ka), bounds, Roots.A42(), rtol=1e-5)
            else
                E[i, j] = E[N-i+2, j]
            end
            # eigenfunctions
            for m = 1:6
                M = system_matrix(E[i, j]; φ, ka, uh)
                uh.c[m, i, j] = (iseven(m) ? 1 : -1) * det(M[2:6, [1:m-1; m+1:6]])
            end
            # normalise coefficients
            κ = [sqrt(E[i, j] - U*cos(φ + 2π*n/3)) for n = 0:2]
            c = view(uh.c, :, i, j)
            X = (c⋅c)a/6 + (c[2]^2 - c[1]^2)sin(2κ[1]a/3)/4κ[1] + (1-cos(2κ[1]a/3))c[1]c[2]/2κ[1] +
                sin(κ[2]a/3)/κ[2] * (c[3]c[4]sin(κ[2]a) + (c[4]^2 - c[3]^2)*cos(κ[2]a)/2) +
                sin(κ[3]a/3)/κ[3] * (c[5]c[6]sin(5κ[3]a/3) + (c[6]^2 - c[5]^2)*cos(5κ[3]a/3)/2)
            c ./= √(N*X) 
        end
    end
end

"""
Construct energy eigenfunctions for each eigenstate in `uh` at each phase number in `whichphases`.
`n_x` specifies the number of points to use for each site.
Return (`x`, `ψ`), where `x` are the abscissas, and `ψ[:, j, i]` = `j`th eigenfunction at `i`th phase.
"""
function make_eigenfunctions(uh::UnperturbedHamiltonian, n_x::Integer, whichphases::AbstractVector{<:Integer})
    (;N, a, U, φₓ, E, c) = uh
    x = range(0, a*N, 3N*n_x+1)
    ψ = Array{ComplexF64,3}(undef, length(x), N, length(whichphases))
    for (i, iϕ) in enumerate(whichphases)
        for j = 1:N # for each eigenenergy
            κ = [sqrt(E[j, iϕ] - U*cos(φₓ[iϕ] + 2π*n/3)) for n = 0:2]
            # construct wave function in the 3 sites of the first cell
            for n = 1:3
                mask = (n-1)*n_x+1:n*n_x
                @. ψ[mask, j, i] = c[2n-1, j, iϕ]sin(κ[n]x[mask]) + c[2n, j, iϕ]cos(κ[n]x[mask])
            end
            # repeat the first cell according to the Bloch's theorem
            for n in 2:N
                @. ψ[3n_x*(n-1)+1:3n_x*n, j, i] = ψ[3n_x*(n-2)+1:3n_x*(n-1), j, i] * cis(2π*j/(N*a))
            end
            ψ[end, j, i] = ψ[1, j, i] # close the loop vor visual convenience
        end
    end
    return x, ψ
end

end