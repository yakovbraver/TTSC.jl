module DeltaModel

import Roots
using LinearAlgebra: eigen, schur, ⋅, svd

"A type for storing the Wannier functions."
mutable struct Wanniers
    E::Matrix{Float64} # `E[j, i]` = mean energy of `j`th wannier at `i`th phase
    pos::Matrix{Float64} # `pos[j, i]` = position eigenvalue of `j`th wannier at `i`th phase
    d::Array{ComplexF64, 3} # `d[:, :, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details
end

"Default-construct an empty `Wanniers` object."
Wanniers() = Wanniers(Float64[;;], Float64[;;], ComplexF64[;;;])

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
    κ::Array{Float64, 3}    # `κ[n, i, j]` = √(E[i, j] - U*cos(φₓ[j] + 2π*n/3))`
    w::Wanniers
end

"""
Construct an `UnperturbedHamiltonian` object.
"""
function UnperturbedHamiltonian(n_cells::Integer; a::Real, λ::Real, U::Real, isperiodic::Bool, φₓ::AbstractVector{<:Real})
    E = Matrix{Float64}(undef, n_cells, length(φₓ))
    c = Array{ComplexF64, 3}(undef, 6, n_cells, length(φₓ))
    κ = Array{Float64, 3}(undef, 3, n_cells, length(φₓ))
    w = Wanniers(Matrix{Float64}(undef, n_cells, length(φₓ)), Matrix{Float64}(undef, n_cells, length(φₓ)),
                 Array{ComplexF64,3}(undef, n_cells, n_cells, length(φₓ)))
    UnperturbedHamiltonian(Int(n_cells), Float64(a), Float64(λ), Float64(U), isperiodic, collect(Float64, φₓ), E, c, κ, w)
end

"Return 𝑔ₙ(𝑥) which realises the pumping protocol."
function 𝑔(x; n, a)
    Int( n/3 <= (x % a)/a < (n+1)/3 )
end

"Return cos(𝑘𝑎) for a given energy `ε` and phase `φ`."
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

"Return the 6×6 matrix that characterises the system."
function system_matrix(uh::UnperturbedHamiltonian, state::Integer, iϕ::Integer)
    (;N, a, λ) = uh
    κ = view(uh.κ, :, state, iϕ)
    s = [sin(n*κ[m]*a/3) for n = 1:3, m = 1:3]
    c = [cos(n*κ[m]*a/3) for n = 1:3, m = 1:3]
    e = cis(-2π*(state-1)/N)
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
                E[i, j] = Roots.find_zero(ε -> cos_ka(ε; φ, uh) - cos(ka), bounds, Roots.A42(), rtol=1e-8)
            else
                E[i, j] = E[N-i+2, j]
            end
            # eigenfunctions
            M = system_matrix(uh, i, j)
            uh.c[:, i, j] = svd(M).V[:, end]
            
            # normalise coefficients
            uh.κ[:, i, j] .= [√(E[i, j] - U*cos(φ + 2π*n/3)) for n = 0:2] # calculate and save κ's
            κ = view(uh.κ, :, i, j)
            c = view(uh.c, :, i, j)
            X = (c⋅c)a/6 + real(c[1]c[2]')sin(κ[1]a/3)^2/κ[1] + (abs2(c[2]) - abs2(c[1]))sin(2κ[1]a/3)/4κ[1] +
                sin(κ[2]a/3)/κ[2] * (real(c[3]c[4]')sin(κ[2]a) + (abs2(c[4]) - abs2(c[3]))cos(κ[2]a)/2) +
                sin(κ[3]a/3)/κ[3] * (real(c[5]c[6]')sin(5κ[3]a/3) + (abs2(c[6]) - abs2(c[5]))cos(5κ[3]a/3)/2)
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
    (;N, a, c, κ) = uh
    x = range(0, a*N, 3N*n_x+1)
    ψ = Array{ComplexF64,3}(undef, length(x), N, length(whichphases))
    for (i, iϕ) in enumerate(whichphases)
        for j = 1:N # for each eigenenergy
            # construct wave function in the 3 sites of the first cell
            for n = 1:3
                mask = (n-1)*n_x+1:n*n_x
                @. ψ[mask, j, i] = c[2n-1, j, iϕ]sin(κ[n, j, iϕ]x[mask]) + c[2n, j, iϕ]cos(κ[n, j, iϕ]x[mask])
            end
            # repeat the first cell according to the Bloch's theorem
            for n in 1:N-1
                @. ψ[3n_x*n+1:3n_x*(n+1), j, i] = ψ[1:3n_x, j, i] * cis(2π*(j-1)n/N)
            end
            ψ[end, j, i] = ψ[1, j, i] # close the loop vor visual convenience
        end
    end
    return x, ψ
end

"Calculate Wannier vectors for the unperturbed Hamiltonian `uh`."
function compute_wanniers!(uh::UnperturbedHamiltonian)
    (;N, a, φₓ, c, E, κ) = uh

    X = Matrix{ComplexF64}(undef, N, N) # position operator
    
    k₂ = 2π/(N*a)
    𝐹(x, i, n, j′, j, iϕ) = begin
        κʲ = κ[i, j, iϕ]
        κʲ′ = κ[i, j′, iϕ]
        cis((n-1)*2π*(j-j′)/N - (κʲ′ + κʲ - k₂)x) / 4 * (
            im * (c[2i-1, j′, iϕ] + im*c[2i, j′, iϕ])' * (c[2i-1, j, iϕ] - im*c[2i, j, iϕ]) / (κʲ′ + κʲ - k₂) +
            (c[2i-1, j, iϕ] + im*c[2i, j, iϕ]) * cis(2κʲ*x) * ( 
                (c[2i, j′, iϕ] - im*c[2i-1, j′, iϕ]) / (-κʲ′ + κʲ + k₂) +
                (c[2i, j′, iϕ] + im*c[2i-1, j′, iϕ]) / ( κʲ′ + κʲ + k₂) * cis(-2κʲ′*x) )' +
            (c[2i-1, j′, iϕ] - im*c[2i, j′, iϕ])' * (c[2i, j, iϕ] + im*c[2i-1, j, iϕ]) * cis(2κʲ′*x) / (κʲ′ - κʲ + k₂) )
    end

    for iϕ in eachindex(φₓ)
        for j in 1:N
            for j′ in 1:N
                X[j′, j] = 0
                for n = 1:N, i = 1:3
                    X[j′, j] += 𝐹((n-1)a + i*a/3, i, n, j′, j, iϕ) - 𝐹((n-1)a + (i-1)a/3, i, n, j′, j, iϕ)
                end
            end
        end
        # `eigen` does not guarantee orthogonality of eigenvectors in case of degeneracies for `X` unitary, so use `schur`
        # (although a degeneracy of coordinates eigenvalues is unlikely here)
        _, uh.w.d[:, :, iϕ], pos_complex = schur(X)
        pos_real = @. (angle(pos_complex) + pi) / k₂ # shift angle from [-π, π) to [0, 2π)
        sp = sortperm(pos_real)                         # sort the eigenvalues
        uh.w.pos[:, iϕ] = pos_real[sp]
        @views Base.permutecols!!(uh.w.d[:, :, iϕ], sp) # sort the eigenvectors in the same way
        uh.w.E[:, iϕ] = [abs2.(dˣ) ⋅ E[:, iϕ] for dˣ in eachcol(uh.w.d[:, :, iϕ])]
    end
end

"""
Construct Wannier functions at each phase number in `whichphases`. `n_x` specifies the number of points to use for each site.
All Wannier functions contained in `uh` are constructed. In the process, energy eigenfunctions are also constructed.
Return `x, ψ, w`, where `ψ[:, j, i]` = `j`th eigenfunction at `i`th phase, and `w[:, j, i]` = `j`th Wannier function at `i`th phase.
"""
function make_wannierfunctions(uh::UnperturbedHamiltonian, n_x::Integer, whichphases::AbstractVector{<:Integer})
    x, ψ = make_eigenfunctions(uh, n_x, whichphases)
    w = similar(ψ)
    for i in eachindex(whichphases)
        for j in 1:uh.N
            w[:, j, i] = sum(uh.w.d[k, j, i] * ψ[:, k, i] for k = 1:uh.N)
        end
    end
    return x, ψ, w
end

end