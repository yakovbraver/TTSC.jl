module DeltaModel

using ProgressMeter: @showprogress
import IntervalRootFinding as RF
using IntervalArithmetic.Symbols # contains `..`
using LinearAlgebra: eigvals, eigvecs, eigen, schur, ⋅, dot, svd, diag, diagm, diagind, Diagonal, Hermitian
using FLoops: @floop, @init
import Optim

"A type for storing the Wannier functions."
mutable struct Wanniers
    targetband::Int
    E::Array{Float64, 3} # `E[j, b, i]` = mean energy of `j`th wannier of the `b`th subband (1 ≤ b ≤ 3) at `i`th phase
    pos::Array{Float64, 3} # `pos[j, b, i]` = position eigenvalue of `j`th wannier of the `b`th subband (1 ≤ b ≤ 3) at `i`th phase
    d::Array{ComplexF64, 4} # `d[:, :, b, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details
end

"""
A type representing the unperturbed Hamiltonian
    ℎ = 𝑝² + 𝜆 ∑ₙ𝛿(𝑥 - 𝑛𝑎/3) + 𝑈 ∑ₙ𝑔ₙ(𝑥)cos(𝜑ₓ + 2π𝑛/3).
"""
mutable struct UnperturbedHamiltonian
    N::Int # number of lattice cells
    a::Float64
    λ::Float64
    U::Float64
    φₓ::Vector{Float64}
    E::Array{Float64, 3}    # `E[ik, b, j]` = `ik`th eigenvalue of `b`th subband at `j`th phase
    c::Array{ComplexF64, 4} # `c[:, ik, b, j]` = `ik`th eigenvector of `b`th subband at `j`th phase
    κ::Array{Float64, 4}    # `κ[n, ik, b, j]` = `√(E[ik, b, j] - U*cos(φₓ[j] + 2π*n/3))`
    w::Wanniers
end

"Construct an `UnperturbedHamiltonian` object."
function UnperturbedHamiltonian(n_cells::Integer; a::Real, λ::Real, U::Real, φₓ::AbstractVector{<:Real})
    E = Float64[;;;]
    c = ComplexF64[;;;;]
    κ = Float64[;;;;]
    w = Wanniers(0, Array{Float64,3}(undef, n_cells, 3, length(φₓ)), Array{Float64,3}(undef, n_cells, 3, length(φₓ)),
                 Array{ComplexF64,4}(undef, n_cells, n_cells, 3, length(φₓ)))
    UnperturbedHamiltonian(Int(n_cells), Float64(a), Float64(λ), Float64(U), collect(Float64, φₓ), E, c, κ, w)
end

"Return 𝑔ₙ(𝑥) which realises the pumping protocol."
function 𝑔(x; n, a)
    Int( n/3 ≤ (x % a)/a < (n+1)/3 )
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

"Same as [`cos_ka`](@ref) calculated using transfer matrix formalism."
function cos_ka_tm(ε::Real; φ::Real, uh::UnperturbedHamiltonian)
    (;a, U, λ) = uh
    κ = [√(ε - U*cos(φ + 2π*n/3)) for n = 0:2]
    R = ComplexF64[1 0; 0 1]
    for n in 0:2
        F = (κ[(n+1)%3 + 1] - im*λ + κ[n+1])cis(κ[n+1]a/3)
        G = (κ[(n+1)%3 + 1] + im*λ - κ[n+1])cis(κ[n+1]a/3)
        R .= [F G'; G F']/2κ[(n+1)%3 + 1] * R
    end 
    return real(R[1, 1])
end

"Return the 6×6 matrix that characterises the system."
function system_matrix(uh::UnperturbedHamiltonian, ik::Integer, sb::Integer, iφ::Integer)
    (;N, a, λ) = uh
    κ = view(uh.κ, :, ik, sb, iφ)
    s = [sin(n*κ[m]*a/3) for n = 1:3, m = 1:3]
    c = [cos(n*κ[m]*a/3) for n = 1:3, m = 1:3]
    e = cis(-2π*(ik-1)/N)
    [0                    -1                  0                    0                   s[3,3]e       c[3,3]e;
     s[1,1]               c[1,1]              -s[1,2]              -c[1,2]             0             0;
     0                    0                   s[2,2]               c[2,2]              -s[2,3]       -c[2,3];
     κ[1]                 -λ                  0                    0                   -κ[3]c[3,3]e  κ[3]s[3,3]e;
     -κ[1]c[1,1]-s[1,1]λ  κ[1]s[1,1]-c[1,1]λ  κ[2]c[1,2]           -κ[2]s[1,2]         0             0;
     0                    0                   -κ[2]c[2,2]-s[2,2]λ  κ[2]s[2,2]-c[2,2]λ  κ[3]c[2,3]    -κ[3]s[2,3]]
end

"""
Diagonalise the unperturbed Hamiltonian `uh`: Find allowed energies at each phase for all bands bracketed in energy by `bounds`,
and calculate the corresponding eigenfunctions.
"""
function diagonalise!(uh::UnperturbedHamiltonian; bounds::Tuple{<:Real, <:Real})
    (;N, a, U, φₓ) = uh

    # perform root-finding once to find the number of subbands in the bounded interval
    rts = RF.roots(ε -> cos_ka(ε; φ=φₓ[1], uh), bounds[1]..bounds[2])
    n_subbands = length(rts)

    uh.E = Array{Float64, 3}(undef, N, n_subbands, length(φₓ)) # the number of different 𝑘's is equal to `N`
    uh.c = Array{ComplexF64, 4}(undef, 6, N, n_subbands, length(φₓ))
    uh.κ = Array{Float64, 4}(undef, 3, N, n_subbands, length(φₓ))

    @showprogress for (iφ, φ) in enumerate(φₓ)
        for ik = 1:N
            ka = 2π*(ik-1)/N

            if ik <= N÷2 + 1
                rts = RF.roots(ε -> cos_ka(ε; φ, uh) - cos(ka), bounds[1]..bounds[2])
                uh.E[ik, :, iφ] = rts .|> RF.root_region .|> RF.mid |> sort
            else
                uh.E[ik, :, iφ] = uh.E[N-ik+2, :, iφ]
            end

            for sb in 1:n_subbands
                # eigenenergies
                uh.κ[:, ik, sb, iφ] .= [√(uh.E[ik, sb, iφ] - U*cos(φ + 2π*n/3)) for n = 0:2]
                
                # eigenfunctions
                M = system_matrix(uh, ik, sb, iφ)
                uh.c[:, ik, sb, iφ] = svd(M).V[:, end]
                
                # normalise coefficients
                κ = view(uh.κ, :, ik, sb, iφ)
                c = view(uh.c, :, ik, sb, iφ)
                X = (c⋅c)a/6 + real(c[1]c[2]')sin(κ[1]a/3)^2/κ[1] + (abs2(c[2]) - abs2(c[1]))sin(2κ[1]a/3)/4κ[1] +
                    sin(κ[2]a/3)/κ[2] * (real(c[3]c[4]')sin(κ[2]a) + (abs2(c[4]) - abs2(c[3]))cos(κ[2]a)/2) +
                    sin(κ[3]a/3)/κ[3] * (real(c[5]c[6]')sin(5κ[3]a/3) + (abs2(c[6]) - abs2(c[5]))cos(5κ[3]a/3)/2)
                c ./= √(N*X) # note division by `N`, which makes the wavefunctions normalised over the whole crystal, not just one cell
            end
        end
    end
end

"""
Construct energy eigenfunctions for each 𝑘 in the `whichband` in `uh`, at each phase number in `whichphases`.
`n_x` specifies the number of points to use for each site.
Return (`x`, `ψ`), where `x` are the abscissas, and `ψ[:, ik, b, i]` = `ik`th eigenfunction of `b`th subband at `i`th phase.
1 ≤ `b` ≤ 3 numbers the subbands of `whichband`.
"""
function make_eigenfunctions(uh::UnperturbedHamiltonian, n_x::Integer, whichband::Integer, whichphases::AbstractVector{<:Integer})
    (;N, a, c, κ) = uh
    x = range(0, a*N, 3N*n_x+1)
    ψ = Array{ComplexF64, 4}(undef, length(x), N, 3, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for b in 1:3
            m = 3(whichband-1) + b
            for ik = 1:N # for each 𝑘
                # construct wave function in the 3 sites of the first cell
                for n = 1:3
                    mask = (n-1)*n_x+1:n*n_x
                    @. ψ[mask, ik, b, i] = c[2n-1, ik, m, iφ]sin(κ[n, ik, m, iφ]x[mask]) + c[2n, ik, m, iφ]cos(κ[n, ik, m, iφ]x[mask])
                end
                # repeat the first cell according to the Bloch's theorem
                for n in 1:N-1
                    @. ψ[3n_x*n+1:3n_x*(n+1), ik, b, i] = ψ[1:3n_x, ik, b, i] * cis(2π*(ik-1)n/N)
                end
                ψ[end, ik, b, i] = ψ[1, ik, b, i] # close the loop vor visual convenience
            end
        end
    end
    return x, ψ
end

"A helper function for calculating ∫𝜓̄ᵢexp(i𝑥𝑘₂)𝜓ⱼ d𝑥."
function 𝐹(uh::UnperturbedHamiltonian, i::Integer, ik′::Integer, ik::Integer, m′::Integer, m::Integer, iφ::Integer, k₂::Real)
    (;κ, c, a) = uh
    κᵐ  = κ[i, ik, m, iφ]
    κᵐ′ = κ[i, ik′, m′, iφ]
    A  = c[2i-1, ik, m, iφ]
    A′ = c[2i-1, ik′, m′, iφ]
    B  = c[2i, ik, m, iφ]
    B′ = c[2i, ik′, m′, iφ]
    F(x) = -cis(-(κᵐ′ + κᵐ - k₂)x)/4 * ( im * (A′ + im*B′)' * (A - im*B) / (κᵐ′ + κᵐ - k₂) +
            (A + im*B) * cis(2κᵐ*x) * ( (B′ - im*A′)/(-κᵐ′ + κᵐ + k₂) + (B′ + im*A′)/(κᵐ′ + κᵐ + k₂)*cis(-2κᵐ′*x) )' +
            (A′ - im*B′)' * (B + im*A) * cis(2κᵐ′*x) / (κᵐ′ - κᵐ + k₂) )
    return F(i*a/3) - F((i-1)*a/3)
end

"A helper function for calculating ∫𝜓̄ᵢ𝜓ⱼ d𝑥."
function 𝐺(uh::UnperturbedHamiltonian, i::Integer, ik′::Integer, ik::Integer, m::Integer, iφ::Integer)
    (;a, κ, c) = uh
    κᵐ = κ[i, ik, m, iφ]
    A  = c[2i-1, ik, m, iφ]
    A′ = c[2i-1, ik′, m, iφ]
    B  = c[2i, ik, m, iφ]
    B′ = c[2i, ik′, m, iφ]
    return 1/2κᵐ * sin(a*κᵐ/3) * ((B′' * B - A′' * A) * cos(a*(2(i-1)+1)*κᵐ/3) +
                                  (A′' * B + B′' * A) * sin(a*(2(i-1)+1)*κᵐ/3)) +
           a/6 * (A′' * A + B′' * B)
       
end

"Calculate Wannier vectors for the unperturbed Hamiltonian `uh`."
function compute_wanniers!(uh::UnperturbedHamiltonian, targetband::Integer)
    (;N, a, φₓ, E) = uh
    uh.w.targetband = targetband

    X = Matrix{ComplexF64}(undef, N, N) # position operator
    
    k₂ = 2π/(N*a)
   
    for iφ in eachindex(φₓ)
        for b in 1:3 # for each of the 3 subbands in the target band
            X .= 0
            m = 3(targetband-1) + b # "global" subband number
            for ik in 1:N
                ik′ = ik % N + 1
                for i = 1:3
                    X[ik′, ik] += 𝐹(uh, i, ik′, ik, m, m, iφ, k₂)
                end
                X[ik′, ik] *= N
            end
            _, uh.w.d[:, :, b, iφ], pos_complex = schur(X)
            pos_real =  @. mod2pi(angle(pos_complex)) / k₂ # shift angle from [-π, π) to [0, 2π)
            sp = sortperm(pos_real)                         # sort the eigenvalues
            uh.w.pos[:, b, iφ] = pos_real[sp]
            @views Base.permutecols!!(uh.w.d[:, :, b, iφ], sp) # sort the eigenvectors in the same way
            uh.w.E[:, b, iφ] = transpose(E[:, m, iφ]) * abs2.(uh.w.d[:, :, b, iφ])
        end
    end
end

"""
Calculate Wannier vectors for the unperturbed Hamiltonian `h` at phase `iφ₀` by mixing the three subbands of the `targetband`.
Return `d, pos, E` as contained in `Wanniers` struct, except that these do not contain a separate dimension for the different subbands.
"""
function compute_wanniers(uh::UnperturbedHamiltonian; iφ₀::Integer=1, targetband::Integer)
    (;N, a, E) = uh

    X = zeros(ComplexF64, 3N, 3N) # position operator
    
    k₂ = 2π/(N*a)

    for b in 1:3  # `b` and `b′` run over the 3 subbands
        m = 3(targetband-1) + b # "global" subband number
        for b′ in 1:3
            m′ = 3(targetband-1) + b′
            for ik in 1:N
                ik′ = ik % N + 1
                for i in 1:3
                    X[N*(b′-1)+ik′, N*(b-1)+ik] += 𝐹(uh, i, ik′, ik, m′, m, iφ₀, k₂)
                end
                X[N*(b′-1)+ik′, N*(b-1)+ik] *= N
            end
        end
    end
    _, d, pos_complex = schur(X) 
    pos_real = @. mod2pi(angle(pos_complex)) / k₂ # shift angle from [-π, π) to [0, 2π)
    sp = sortperm(pos_real)          # sort the eigenvalues
    pos = pos_real[sp]
    @views Base.permutecols!!(d, sp) # sort the eigenvectors in the same way
    E = transpose(uh.E[range((iφ₀-1)size(uh.E, 2)size(uh.E, 1) + 3(targetband-1)size(uh.E, 1) + 1, length=3N)]) * abs2.(d)
    return d, pos, vec(E) # materialise `E`
end

"""
Construct Wannier functions at each phase number in `whichphases`. `n_x` specifies the number of points to use for each site.
All Wannier functions contained in `uh` are constructed. In the process, energy eigenfunctions are also constructed.
Return `x, ψ, w`, where `ψ[:, ik, b, i]` = `ik`th eigenfunction of `b`the subband at `i`th phase,
and `w[:, j, b, i]` = `j`th Wannier function of `b`th subband at `i`th phase.
"""
function make_wannierfunctions(uh::UnperturbedHamiltonian, n_x::Integer, whichphases::AbstractVector{<:Integer})
    x, ψ = make_eigenfunctions(uh, n_x, uh.w.targetband, whichphases)
    w = similar(ψ)
    for (i, iφ) in enumerate(whichphases)
        for b in 1:3
            w[:, :, b, i] = ψ[:, :, b, i] * uh.w.d[:, :, b, iφ]
        end
    end
    return x, ψ, w
end

"Abstract supertype of the tight-binding Hamiltonians."
abstract type AbstractTBHamiltonian end

"""
A type representing the tight-binding Hamiltonian
    ℎₜ(𝜑ₓ) = ∑ⱼ 𝐽₁𝑏⁺ⱼ𝑎ⱼ + 𝐽₂𝑐⁺ⱼ𝑏ⱼ + 𝐽₃𝑎⁺ⱼ₊₁𝑐ⱼ + h.c.
             + 𝑈 ∑ⱼ 𝑎⁺ⱼ𝑎ⱼcos(𝜑ₓ) + 𝑏⁺ⱼ𝑏ⱼcos(𝜑ₓ + 2π/3) + 𝑐⁺ⱼ𝑐ⱼcos(𝜑ₓ + 4π/3)
"""
mutable struct SimpleTBHamiltonian <: AbstractTBHamiltonian
    N::Int # number of lattice cells
    a::Float64
    U::Float64
    J::Vector{ComplexF64}
    isperiodic::Bool
    φₓ::Vector{Float64}
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, 1 ≤ i ≤ 3N, 1 ≤ j ≤ length(φₓ)
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::Wanniers 
end

"Construct a `SimpleTBHamiltonian` object."
function SimpleTBHamiltonian(n_cells::Integer; a::Real, U::Real, J::Vector{<:Number}, isperiodic::Bool, φₓ::AbstractVector{<:Real})
    E = Matrix{Float64}(undef, 3n_cells, length(φₓ))
    c = Array{ComplexF64, 3}(undef, 3n_cells, 3n_cells, length(φₓ))
    w = Wanniers(0, Array{Float64,3}(undef, n_cells, 3, length(φₓ)), Array{Float64,3}(undef, n_cells, 3, length(φₓ)),
                 Array{ComplexF64,4}(undef, n_cells, n_cells, 3, length(φₓ)))
    SimpleTBHamiltonian(Int(n_cells), Float64(a), Float64(U), ComplexF64.(J), isperiodic, collect(Float64, φₓ), E, c, w)
end

"Diagonalise the TB Hamiltonian `tbh` at each phase."
function diagonalise!(tbh::SimpleTBHamiltonian)
    (;N, U, J, isperiodic, φₓ) = tbh
    for (i, φ) in enumerate(φₓ)
        diag = repeat([U*cos(φ), U*cos(φ+2π/3), U*cos(φ+4π/3)], N)
        J_diag = [repeat(J, N-1); J[1:2]]
        H = diagm(0 => diag, -1 => J_diag, 1 => conj.(J_diag))
        if isperiodic
            H[1, end] = J[3]
            H[end, 1] = J[3]'
        end
        tbh.E[:, i], tbh.c[:, :, i] = eigen(Hermitian(H))
    end
end

"Return the 𝑘-space Hamiltonian matrix for `tbh` at the given phase `φ` and at 𝑘𝑎 = `ka`."
function kspace_hamiltonian(tbh::SimpleTBHamiltonian, φ::Real, ka::Real)
    (;U, J) = tbh
    [U*cos(φ)        J[1]'         J[3]cis(-ka)
     J[1]            U*cos(φ+2π/3) J[2]'
     (J[3]cis(-ka))' J[2]          U*cos(φ+4π/3)]
end

"""
Diagonalise the TB Hamiltonian `tbh` in 𝑘-space at each phase for the values of 𝑘𝑎 in `ka`.
Return the matrix of eigenenergies `E`, where `E[:, i]` is the energy at `i`th phase.
In `E`, rows 1:3 corresopnd to `ka[1]`, rows 4:6 correspond to `ka[2]`, and so on.
"""
function diagonalise_kspace(tbh::SimpleTBHamiltonian, ka::AbstractVector{<:Real})
    E = Matrix{Float64}(undef, 3length(ka), length(tbh.φₓ))
    for (iφ, φ) in enumerate(tbh.φₓ)
        for ik in eachindex(ka)
            E[3(ik-1)+1:3ik, iφ] .= eigvals(kspace_hamiltonian(tbh, φ, ka[ik]))
        end
    end
    return E
end

""" 
A type representing a general tight-binding Hamiltonian
    ℎₜ(𝜑ₓ) = ∑ᵢⱼ 𝐽ᵢⱼ(𝜑ₓ)𝑎⁺ᵢ𝑎ⱼ
"""
mutable struct TBHamiltonian <: AbstractTBHamiltonian
    N::Int # number of lattice cells
    a::Float64
    H::Array{ComplexF64, 3} # Hamiltonian matrix
    isperiodic::Bool
    φₓ::Vector{Float64}
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, 1 ≤ i ≤ 3N, 1 ≤ j ≤ length(φₓ)
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::Wanniers 
end

"Construct a `TBHamiltonian` for the `targetband` of unperturbed Hamiltonian `uh`. `d` is a matrix of Wannier vectors construced at phase `iφ₀`."
function TBHamiltonian(uh::UnperturbedHamiltonian; d::Matrix{ComplexF64}, iφ₀::Integer=1, isperiodic::Bool, targetband::Integer)
    (;a, N, U, φₓ) = uh
    n_φₓ = length(φₓ)
    n_w = size(d, 1) # number of Wanniers
    H = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ) # TB Hamiltonian matrix

    H₀ = d' * Diagonal(uh.E[range((iφ₀-1)size(uh.E, 2)size(uh.E, 1) + N*3(targetband-1) + 1, length=n_w)]) * d

    ψ∑ψ = Matrix{ComplexF64}(undef, n_w, n_w)

    for (iφ, φ) in enumerate(φₓ)
        ψ∑ψ .= 0
        for b in 1:3  # `b` and `b′` run over the 3 subbands
            m = 3(targetband-1) + b # "global" subband number
            for b′ in 1:3
                m′ = 3(targetband-1) + b′
                for ik in 1:N
                    if m′ == m
                        ψ∑ψ[N*(b′-1)+ik, N*(b-1)+ik] = sum(𝐺(uh, r, ik, ik, m, iφ₀) * (cos(φ + 2π*(r-1)/3) - cos(φₓ[iφ₀] + 2π*(r-1)/3)) for r in 1:3)
                    else
                        for r in 1:3
                            ψ∑ψ[N*(b′-1)+ik, N*(b-1)+ik] += 𝐹(uh, r, ik, ik, m′, m, iφ₀, 0) * (cos(φ + 2π*(r-1)/3) - cos(φₓ[iφ₀] + 2π*(r-1)/3))
                        end
                    end
                    ψ∑ψ[N*(b′-1)+ik, N*(b-1)+ik] *= N
                end
            end
        end
        H[:, :, iφ] = H₀ + U * d' * ψ∑ψ * d
    end

    !isperiodic && (H[1, end, :] .= H[end, 1, :] .= 0)

    E = Matrix{Float64}(undef, n_w, n_φₓ)
    c = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ)
    w = Wanniers(0, Array{Float64,3}(undef, N, 3, length(φₓ)), Array{Float64,3}(undef, N, 3, length(φₓ)),
                 Array{ComplexF64,4}(undef, N, N, 3, length(φₓ)))
    TBHamiltonian(N, a, H, isperiodic, φₓ, E, c, w)
end

"Diagonalise the TB Hamiltonian `tbh` at each phase."
function diagonalise!(tbh::TBHamiltonian)
    for iφ in eachindex(tbh.φₓ)
        tbh.E[:, iφ], tbh.c[:, :, iφ] = eigen(Hermitian(tbh.H[:, :, iφ]))
    end
end

"Calculate Wannier vectors for each of the three subbands for the TB Hamiltonian `tbh`."
function compute_wanniers!(tbh::AbstractTBHamiltonian)
    (;N, φₓ) = tbh
    for b in 1:3
        levels = N*(b-1)+1:N*b
        if tbh.isperiodic
                X = Diagonal(cis.(2π/3N .* (0:3N-1))) # position operator in coordinate representation
            for iφ in eachindex(φₓ)
                XE = tbh.c[:, levels, iφ]' * X * tbh.c[:, levels, iφ] # position operator in energy representation
                _, tbh.w.d[:, :, b, iφ], pos_complex = schur(XE)
                pos_real = @. mod2pi(angle(pos_complex)) / 2π * N # shift angle from [-π, π) to [0, 2π)
                sp = sortperm(pos_real)                        # sort the eigenvalues
                tbh.w.pos[:, b, iφ] = pos_real[sp]
                @views Base.permutecols!!(tbh.w.d[:, :, b, iφ], sp) # sort the eigenvectors in the same way
                tbh.w.E[:, b, iφ] = transpose(tbh.E[levels, iφ]) * abs2.(tbh.w.d[:, :, b, iφ])
            end
        else
            X = Diagonal((0:3N-1)./3) # position operator in coordinate representation
            for iφ in eachindex(φₓ)
                XE = tbh.c[:, levels, iφ]' * X * tbh.c[:, levels, iφ] # position operator in energy representation
                tbh.w.pos[:, b, iφ], tbh.w.d[:, :, b, iφ] = eigen(Hermitian(XE))
                tbh.w.E[:, b, iφ] = transpose(tbh.E[levels, iφ]) * abs2.(tbh.w.d[:, :, b, iφ])
            end
        end
    end
end

"""
Construct Wannier functions at each phase number in `whichphases`. All Wannier functions contained in `tbh` are constructed.
Return `w`, where `w[:, j, b i]` = `j`th Wannier function of `b`th subband at `i`th phase.
"""
function make_wannierfunctions(tbh::AbstractTBHamiltonian, whichphases::AbstractVector{<:Integer})
    (;N) = tbh
    w = Array{ComplexF64, 4}(undef, size(tbh.c, 1), N, 3, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for b in 1:3
            w[:, :, b, i] = tbh.c[:, range(N*(b-1)+1, length=N), iφ] * tbh.w.d[:, :, b, iφ]
        end
    end
    return w
end

"A type for storing the Floquet Wannier functions."
mutable struct FloquetWanniers
    targetsubbands::Vector{Int} # in the TB case, these are actually target levels, not subbands
    E::Array{Float64, 2} # `E[j, i]` = mean energy of `j`th wannier at `i`th phase
    pos::Array{Float64, 2} # `pos[j, i]` = position eigenvalue of `j`th wannier at `i`th phase
    d::Array{ComplexF64, 3} # `d[:, :, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details
end

"Default-construct an empty `FloquetWanniers` object."
FloquetWanniers() = FloquetWanniers(Int[], Float64[;;], Float64[;;], ComplexF64[;;;])

"Swap energies, positions, and vectors of wanniers `i` and `j` at phases `iφ`."
function swap_wanniers!(w::FloquetWanniers, i::Integer, j::Integer, iφ::AbstractVector{<:Integer})
    temp_E = w.E[i, iφ]
    w.E[i, iφ] = w.E[j, iφ]
    w.E[j, iφ] = temp_E
    
    temp_pos = w.pos[i, iφ]
    w.pos[i, iφ] = w.pos[j, iφ]
    w.pos[j, iφ] = temp_pos
    
    temp_d = w.d[:, i, iφ]
    w.d[:, i, iφ] = w.d[:, j, iφ]
    w.d[:, j, iφ] = temp_d

    return nothing
end

"""
A type representing the Floquet Hamiltonian
    ℋ = ℎ - i∂ₜ + λₛcos(12π𝑥/𝑎)cos(2𝜔𝑡) + λₗcos(6π𝑥/𝑎)cos(𝜔𝑡 + 𝜑ₜ),
where ℎ is the unperturbed Hamiltonian represented by [`UnperturbedHamiltonian`](@ref).
"""
mutable struct FloquetHamiltonian
    uh::UnperturbedHamiltonian
    s::Int
    λₛ::Float64
    λₗ::Float64
    ω::Float64
    pumptype::Symbol
    E::Array{Float64, 3}    # `E[i, ik, j]` = Floquet quasienergy of `i`th subband and `ik`th quasimomentum at `j`th phase
    b::Array{ComplexF64, 4} # `b[:, i, ik, j]` = `i`th eigenvector `i`th subband and `ik`th quasimomentum at `j`th phase
    ν::Vector{Int}  # band map 𝜈(𝑚)
    w::FloquetWanniers
end

"""
Construct a `FloquetHamiltonian` object. `minband` is the first energy band of `uh` to use when constructing the Floquet Hamiltonian matrix.
Type of pumping is controlled via `pumptype`: `:time` for temporal, `:space` for spatial, or anything else for simultaneous space-time pumping.
In the case of time-only pumping, it is assumed that 𝜑ₓ = 0, and hence that `uh.φₓ[1] == 0`.
"""
function FloquetHamiltonian(uh::UnperturbedHamiltonian; s::Integer, λₛ::Real, λₗ::Real, ω::Real, pumptype::Symbol)
    n_levels = size(uh.E, 2)
    ν = [ceil(Int, m/3) for m in 1:n_levels]
    
    E = Array{Float64, 3}(undef, n_levels, uh.N, length(uh.φₓ))
    b = Array{ComplexF64, 4}(undef, n_levels, n_levels, uh.N, length(uh.φₓ))
    
    FloquetHamiltonian(uh, Int(s), Float64(λₛ), Float64(λₗ), Float64(ω), pumptype, E, b, ν, FloquetWanniers())
end

"""
Diagonalise the Floquet Hamiltonian `fh` at each phase. If `reorder=true`, then reorder the eigenstates according to the unperturbed order.
By default, the pumping protocol 𝜑ₜ = 𝜑ₓ is assumed. Alternatively, you can pass a vector `φₜ` of temporal phases, whose length matches that of `fh.uh.φₓ`.
"""
function diagonalise!(fh::FloquetHamiltonian; reorder::Bool=false, φₜ::AbstractVector{<:Real}=fh.uh.φₓ)
    length(fh.uh.φₓ) != length(φₜ) && (@error "length(fh.uh.φₓ) != length(φₜ). Aborting."; return)
    (;N, a, φₓ, E) = fh.uh
    (;s, ω, λₛ, λₗ, pumptype, ν) = fh

    n_levels = size(fh.E, 1)

    H = zeros(ComplexF64, n_levels, n_levels) # ℋ matrix, will only fill the lower triangle

    for ik in 1:N
        for iφ in eachindex(φₓ)
            for m in 1:n_levels
                # for time-only pumping, always take the eigenenergies at the first phase, which is asssumed to correspond to 𝜑ₓ = 0
                p = (pumptype == :time ? 1 : iφ)
                H[m, m] = E[ik, m, p] - ν[m]*ω/s

                # place the elements of the long lattice
                for g in 1:3
                    m′ = 3(s + ν[m] - 1) + g
                    m′ > n_levels && break
                    if pumptype != :time || iφ == 1 # if pumping is time-only, this must be calculated only once, at `iφ` = 1
                        ∫cos = ComplexF64(0)
                        for i = 1:3, k₂ in (-6π/a, 6π/a)
                            ∫cos += 𝐹(fh.uh, i, ik, ik, m′, m, iφ, k₂)
                        end
                        ∫cos *= N # restore proper normalisation; `fh.uh.c` used in `𝐹` are normalised over all the cells, but we need one-cell normalisation here
                        # if pumping is space-time, then also multiply by cis(-𝜑ₜ)
                        H[m′, m] = (pumptype == :space ? λₗ/4 * ∫cos : λₗ/4 * ∫cos * cis(-φₜ[iφ]))
                    elseif pumptype == :time 
                        H[m′, m] *= cis(-(φₜ[iφ]-φₜ[iφ-1]))
                    end
                end
                
                # place the elements of the short lattice
                if pumptype != :time || iφ == 1 # if pumping is time-only, this must be calculated only once, at `iφ` = 1
                    for g in 1:3
                        m′ = 3(2s + ν[m] - 1) + g
                        m′ > n_levels && break
                        ∫cos = ComplexF64(0)
                        for i = 1:3, k₂ in (-12π/a, 12π/a)
                            ∫cos += 𝐹(fh.uh, i, ik, ik, m′, m, iφ, k₂)
                        end
                        ∫cos *= N # restore proper normalisation; `fh.uh.c` used in `𝐹` are normalised over all the cells, but we need one-cell normalisation here
                        H[m′, m] = λₛ/4 * ∫cos
                    end
                end
            end
            fh.E[:, ik, iφ], fh.b[:, :, ik, iφ] = eigen(Hermitian(H, :L))
            if reorder
                perm = diag(H) .|> real |> sortperm |> invperm # get a permutation "restoring" the order of spatial bands
                @views permute!(fh.E[:, ik, iφ], perm)
                @views Base.permutecols!!(fh.b[:, :, ik, iφ], perm)
            end
        end
    end
end

"Swap energies and eigenvectors of states `i` and `j` for quasimomenta `ik` at phases `iφ`."
function swap_eigenstates!(fh::FloquetHamiltonian, i::Integer, j::Integer, ik::AbstractVector{<:Integer}, iφ::AbstractVector{<:Integer})
    temp_E = fh.E[i, ik, iφ]
    fh.E[i, ik, iφ] = fh.E[j, ik, iφ]
    fh.E[j, ik, iφ] = temp_E

    temp_b = fh.b[:, i, ik, iφ]
    fh.b[:, i, ik, iφ] = fh.b[:, j, ik, iφ]
    fh.b[:, j, ik, iφ] = temp_b

    return nothing
end

"""
Construct Floquet modes at coordinates in `x` and time moments in `Ωt` for each spatial subband number in `whichsubband` at each phase number in `whichphases`.
Return `x, u`, where `x` are the abscissas and `u[ix, it, ik, j, i]` = wavefunction corresponding to the `ik`th value of 𝑘 and `j`th spatial subband
at `whichphases[i]`th phase at `ix`th coordinate at `it`th time moment.
"""
function make_eigenfunctions(fh::FloquetHamiltonian, n_x::Integer, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer},
                             whichsubbands::AbstractVector{<:Integer})
    (;N, a) = fh.uh
    x = range(0, a*N, 3N*n_x+1)
    u = Array{ComplexF64, 5}(undef, 3N*n_x+1, length(Ωt), N, length(whichsubbands), length(whichphases))
    n_levels = size(fh.E, 1) # number of Floquet levels; equivalently, number of levels of ℎ
    # Eigenfunctions of ℎ, which are mixed during construction of `u`. For time-only pumping use only eigenstates at the first phase, corresponding to 𝜑ₓ = 0
    ψ = Array{ComplexF64, 5}(undef, 3N*n_x+1, N, 3, n_levels÷3, (fh.pumptype == :time ? 1 : length(whichphases))) # ψ[ix, ik, subband, band, iφ]
    for band in 1:n_levels÷3
        _, ψ[:, :, :, band, :] = make_eigenfunctions(fh.uh, n_x, band, (fh.pumptype == :time ? [1] : whichphases))
    end
    for (i, iφ) in enumerate(whichphases)
        p = (fh.pumptype == :time ? 1 : i)
        for (j, js) in enumerate(whichsubbands)
            for ik in 1:N
                for (it, t) in enumerate(Ωt)
                    u[:, it, ik, j, i] = sum(cis(-fh.ν[m]*t) * ψ[:, ik, (m-1)%3+1, fh.ν[m], p] * fh.b[m, js, ik, iφ] for m in 1:n_levels)
                end
            end
        end
    end
    return x, u
end

"""
Calculate Wannier vectors for the Floquet Hamiltonian `fh` using the quasienergy levels `targetsubbands`.
`Ωt` is the time moment at which to diagonalise the coordinate operator; pass `slide_time=true` if you wish this moment
to change linearly with phase, sliding by π over the period. This may be beneficial for temporal pumping.
"""
function compute_wanniers!(fh::FloquetHamiltonian; targetsubbands::AbstractVector{<:Integer}, Ωt::Real=π/5, slide_time=false)
    (;N, a, φₓ) = fh.uh

    n_w = length(targetsubbands) * N
    E = Matrix{Float64}(undef, n_w, length(φₓ))
    pos = Matrix{Float64}(undef, n_w, length(φₓ))
    d = Array{ComplexF64, 3}(undef, n_w, n_w, length(φₓ))
    fh.w = FloquetWanniers(targetsubbands, E, pos, d)
    
    n_levels = size(fh.E, 1)
    n_subbands = length(targetsubbands)    
    
    k₂ = 2π/(N*a)
    
    @floop for iφ in eachindex(φₓ)
        @init begin
            expik = Array{ComplexF64, 4}(undef, n_levels, n_levels, N, N)
            X = Matrix{ComplexF64}(undef, n_w, n_w) # position operator
            first_iter = true # idicates if the current thread is on its first iteration
        end
        # if pumping is time-only, then `expik` can be calculated only once, using `c`'s at 𝜑ₓ = 0
        if fh.pumptype != :time || first_iter
            expik .= 0
            for ik in 1:N
                ik′ = ik % N + 1
                for m in 1:n_levels,  m′ in 1:n_levels
                    for i in 1:3
                        expik[m′, m, ik′, ik] += 𝐹(fh.uh, i, ik′, ik, m′, m, (fh.pumptype == :time ? 1 : iφ), k₂)
                    end
                    expik[m′, m, ik′, ik] *= N
                end
            end
            first_iter = false
        end

        t = slide_time ? Ωt - (iφ-1)/length(φₓ) * π : Ωt

        for ik in 1:N,  ik′ in 1:N
            for (in, n) in enumerate(targetsubbands)
                for (in′, n′) in enumerate(targetsubbands)
                    X[in′+(ik′-1)*n_subbands, in+(ik-1)*n_subbands] = sum(fh.b[m, n, ik, iφ] * sum(fh.b[m′, n′, ik′, iφ]' * expik[m′, m, ik′, ik] *
                                                                          cis((fh.ν[m′] - fh.ν[m]) * t) for m′ in 1:n_levels) for m in 1:n_levels)
                end
            end
        end
        _, d[:, :, iφ], pos_complex = schur(X)
        pos_real = @. mod2pi(angle(pos_complex)) / k₂ # shift angle from [-π, π) to [0, 2π)
        sp = sortperm(pos_real)         # sort the eigenvalues
        pos[:, iφ] = pos_real[sp]
        @views Base.permutecols!!(d[:, :, iφ], sp) # sort the eigenvectors in the same way
        E[:, iφ] = [sum(abs2(dˣ[m]) * fh.E[targetsubbands[(m-1)%n_subbands+1], (m-1)÷n_subbands+1, iφ] for m in eachindex(dˣ)) for dˣ in eachcol(d[:, :, iφ])]
    end
end

"""
Construct Wannier functions at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `fh` are constructed.
In the process, quasienergy eigenfunctions (Floquet modes) are also constructed.
Return `x, u, w`, where `w[ix, it, j, i]` = `j`th Wannier function at `whichphases[i]`th phase at `ix`th coordinate at `it`th time moment,
`u` is an array of Floquet modes (`u[ix, it, ik, j, i]`), and `x` contains the abscissas.
"""
function make_wannierfunctions(fh::FloquetHamiltonian, n_x::Integer, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
    (;N) = fh.uh
    n_subbands = length(fh.w.targetsubbands)
    n_w = n_subbands * N
    x, u = make_eigenfunctions(fh, n_x, Ωt, whichphases, fh.w.targetsubbands) # format: `u[ix, it, ik, j, i]`
    w = Array{ComplexF64, 4}(undef, length(x), length(Ωt), n_w, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for j in 1:n_w
            w[:, :, j, i] = sum(fh.w.d[m, j, iφ] * u[:, :, (m-1)÷n_subbands+1, (m-1)%n_subbands+1, i] for m = 1:n_w)
        end
    end
    return x, u, w
end

"""
Find the optimal superposition for each pair of Wannier numbers in `whichstates` at phase `iφ`,
such that the overlap of probability density is minimised. Overwrite the corresponding states `fh.w.d`.
Note that `fh.w.E` and `fh.w.pos` are not recalculated.
"""
function optimise_wanniers!(fh::FloquetHamiltonian; whichstates::Vector{<:Tuple{Integer, Integer}}, iφ::Integer)
    n_x = 10
    Ωt = range(0, 2π, length=10fh.s)
    x, _, w = DeltaModel.make_wannierfunctions(fh, n_x, Ωt, [iφ])

    for (i, j) in whichstates
        wᵢ = @view w[:, :, i, 1]
        wⱼ = @view w[:, :, j, 1]
        
        best_x = Vector{Float64}(undef, 3)
        best_val = 0.0
        for _ in 1:10
            opt = Optim.optimize([π * rand(), 2π * rand(), 2π * rand()], Optim.NelderMead()) do x
                local U = get_U(x...)
                s = 0.0
                for k in eachindex(wᵢ)
                    wᵢ′ = U[1] * wᵢ[k] + U[2] * wⱼ[k]
                    wⱼ′ = U[3] * wᵢ[k] + U[4] * wⱼ[k]
                    s -= abs2(abs2(wᵢ′) - abs2(wⱼ′))
                end
                s
            end
            val = Optim.minimum(opt)
            if val < best_val
                best_val = val
                best_x .= Optim.minimizer(opt) 
            end
        end
        U = get_U(best_x...)

        for k in axes(fh.w.d, 1)
            dᵢ = fh.w.d[k, i, iφ]
            fh.w.d[k, i, iφ] = U[1] * dᵢ + U[2] * fh.w.d[k, j, iφ]
            fh.w.d[k, j, iφ] = U[3] * dᵢ + U[4] * fh.w.d[k, j, iφ]
        end
    end
end

"Generate a 2×2 unitary matrix, returned as a tuple, parameterised by the angles 0 ≤ `ϑ` ≤ π, 0 ≤ `ϕ`, `α` ≤ 2π."
function get_U(ϑ::Real, ϕ::Real, α::Real)
    n1, n2, n3 = sin(ϑ)cos(ϕ), sin(ϑ)sin(ϕ), cos(ϑ)
    s, c = sincos(α)
    return c + im * n3 * s,        im * n1 * s + n2 * s,
           im * n1 * s - n2 * s,   c - im * n3 * s
end

"""
For each phase, reorder the Wanniers in each spatial site so that the first Wannier has its "right" probability density maximum in the region
3π/2 < 𝛺𝑡 < 2π, the second has its right maximum in π < 𝛺𝑡 < 3π/2, third in π/2 < 𝛺𝑡 < π, and fourth in 0 < 𝛺𝑡 < π/2.
In practice, this requires either swapping first and second Wanniers, or second with third and then the resulting second with first.
Optionally, run an optimisation of the second and third Wanniers in each spatial site before reordering.
"""
function order_wanniers!(fh::FloquetHamiltonian; optimise::Bool=true)
    (;N, φₓ) = fh.uh
    n_x = 20
    n_t = 20fh.s
    Ωt = range(0, 2π, length=n_t)

    n_w = size(fh.w.d, 2)

    if optimise
        whichstates = [(4(s-1)+2, 4(s-1)+3) for s in 1:n_w÷4]
        @floop for iφ in eachindex(φₓ)
            optimise_wanniers!(fh; whichstates, iφ)
        end
    end
    
    I = Vector{Float64}(undef, 3) # stores the value of ∫|𝑤ᵢ|²d𝑥 d𝑡 over 3𝑙/4 < 𝑥 < 𝑙, 3π/2 < 𝛺𝑡 < 2π for 𝑖 = 1, 2, 3
    _, _, w = make_wannierfunctions(fh, n_x, Ωt, eachindex(φₓ))
    for iφ in eachindex(φₓ)
        for s in 1:n_w÷4 # for each group of 4 Wanniers (i.e. for each spatial site)
            for j in 1:3 # for each of the first 3 Wanniers in a given spatial site
                I[j] = sum(abs2, @view w[(s-1)n_x+3n_x÷4:s*n_x, 7n_t÷10:end, 4(s-1)+j, iφ])
            end
            m = argmax(I)
            if m == 2
                swap_wanniers!(fh.w, 4(s-1)+1, 4(s-1)+2, [iφ])
            elseif m == 3
                swap_wanniers!(fh.w, 4(s-1)+2, 4(s-1)+3, [iφ])
                swap_wanniers!(fh.w, 4(s-1)+1, 4(s-1)+2, [iφ])
            end
        end
    end

    # additional reordering in every other spatial site to get separable TB Hamiltonian
    for i in 2:2:3N
        j = 4(i-1)
        swap_wanniers!(fh.w, j+3, j+1, eachindex(φₓ))
        swap_wanniers!(fh.w, j+4, j+2, eachindex(φₓ))
    end
end

"""
A type representing a 2D (time+space) tight-binding Hamiltonian.
"""
mutable struct TBFloquetHamiltonian
    N::Int
    a::Float64
    H::Array{ComplexF64, 3} # Hamiltonian matrix
    periodicity::Symbol
    φₓ::Vector{Float64}
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, `1 ≤ i ≤ 6N`, `1 ≤ j ≤ length(φₓ)`
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::FloquetWanniers 
end

"""
Experimental: Construct a `TBFloquetHamiltonian` object using the periodic Wanniers at a single phase contained in `fh`.
`fh` should be calculated for a single spatial site. `pumptype` may or may not coincide with `fh.pumptype`.
`startsubband` is the first subband of `fh` to use for constructing Wanniers.
"""
function TBFloquetHamiltonian(fh::FloquetHamiltonian, startsubband::Integer, pumptype::Symbol)
    (;a, U, φₓ) = fh.uh
    (;s, λₗ, ν) = fh
    n_φₓ = length(φₓ)
    n_w = 4 * 3 # number of temporal sites * number of spatial sites
    H = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ)

    n_m = size(fh.E, 1) # number of levels of ℎ considered
    Ψ = Matrix{ComplexF64}(undef, n_m, n_m)

    iφ₀ = 1 # phase index at which to take the Wanniers -- any choice should work, but it is simpler to assume `iφ₀ = 1` below
    ik = 1

    d = @view fh.w.d[:, :, iφ₀]
    bd = fh.b[:, range(startsubband, length=n_w), ik, iφ₀] * d

    H₀ = d' * Diagonal(fh.E[range(startsubband, length=n_w), ik, iφ₀]) * d

    for (iφ, φ) in enumerate(φₓ)
        for m in 1:n_m # `m` is the subband index of ℎ
            for m′ in 1:n_m
                Ψ[m′, m] = 0
                if pumptype != :time # if pumping is not time-only, account for the change of the spatial phase
                    if m′ == m
                        Ψ[m′, m] += U * sum(𝐺(fh.uh, r, ik, ik, m, iφ₀) * (cos(φ + 2π*(r-1)/3) - cos(2π*(r-1)/3)) for r in 1:3)
                    elseif ν[m] == ν[m′]
                        for r in 1:3
                            Ψ[m′, m] += U * 𝐹(fh.uh, r, ik, ik, m′, m, iφ₀, 0) * (cos(φ + 2π*(r-1)/3) - cos(2π*(r-1)/3))
                        end
                    end
                end
                if pumptype != :space # if pumping is not space-only, account for the change of the temporal phase
                    if ν[m] == ν[m′] + s
                        e = cis(+φ)
                    elseif ν[m] == ν[m′] - s
                        e = cis(-φ)
                    else
                        continue
                    end
                    for r = 1:3, k₂ in (-6π/a, 6π/a)
                        Ψ[m′, m] += λₗ/4 * (e - 1) * 𝐹(fh.uh, r, ik, ik, m′, m, iφ₀, k₂)
                    end
                end
            end
        end
        H[:, :, iφ] = H₀ + bd' * Ψ * bd
    end

    E = Matrix{Float64}(undef, n_w, n_φₓ)
    c = Array{ComplexF64, 3}(undef, n_w, n_w, n_φₓ)
    TBFloquetHamiltonian(1, a, H, :both, φₓ, E, c, FloquetWanniers())
end

"""
Construct a `TBFloquetHamiltonian` object using the corresponding Wanniers contained in `fh` at each phase.
`perdicity` should be `:none`, `:space`, `:time`, or anything else for toroidal time-space periodicity.
"""
function TBFloquetHamiltonian(fh::FloquetHamiltonian; periodicity::Symbol)
    (;N, φₓ) = fh.uh
    n_φₓ = length(φₓ)
    n_s = length(fh.w.targetsubbands) # number of spatial subbands of ℋ mixed
    n_t = 4 # number of temporal sites
    H = Array{ComplexF64, 3}(undef, n_s*N, n_s*N, n_φₓ)
    for iφ in eachindex(φₓ)
        for a = 1:n_s*N, b = a:n_s*N
            H[b, a, iφ] = sum(fh.w.d[i, a, iφ]' * fh.w.d[i, b, iφ] * fh.E[fh.w.targetsubbands[(i-1)%n_s+1], (i-1)÷n_s+1, iφ] for i = axes(fh.w.d, 1))
            H[a, b, iφ] = H[b, a, iφ]'
        end
        if periodicity == :space || periodicity == :none # make non-periodic in time
            for i in 0:3N-1
                H[n_t*i+1, n_t*(i+1), iφ] = H[n_t*(i+1), n_t*i+1, iφ] = 0
            end
        end
        if periodicity == :time || periodicity == :none # make non-periodic in space
            H[1:n_t, end-n_t+1:end, iφ] .= H[end-n_t+1:end, 1:n_t, iφ] .= 0
        end
    end

    E = Matrix{Float64}(undef, n_s*N, n_φₓ)
    c = Array{ComplexF64, 3}(undef, n_s*N, n_s*N, n_φₓ)
    TBFloquetHamiltonian(N, fh.uh.a, H, periodicity, φₓ, E, c, FloquetWanniers())
end

"Diagonalise the TB Floquet Hamiltonian `tbh` at each phase."
function diagonalise!(tbh::TBFloquetHamiltonian)
    for i in eachindex(tbh.φₓ)
        tbh.E[:, i], tbh.c[:, :, i] = eigen(Hermitian(tbh.H[:, :, i]))
    end
end

"""
Calculate Wannier vectors for the TB floquet Hamiltonian `tbh` using the energy levels `targetlevels`.
These levels are counted from the spectrum of `tbh`, starting from the lowest.
"""
function compute_wanniers!(tbh::TBFloquetHamiltonian, targetlevels::AbstractVector{<:Integer})
    (;N, φₓ, periodicity) = tbh

    n_w = length(targetlevels)
    n_t = 4 # number of temporal sites
    E = Matrix{Float64}(undef, n_w, length(φₓ))
    pos = Matrix{Float64}(undef, n_w, length(φₓ))
    d = Array{ComplexF64, 3}(undef, n_w, n_w, length(φₓ))
    tbh.w = FloquetWanniers(targetlevels, E, pos, d) # use `targetlevels` in place of `targetsubbands` argument

    if periodicity == :none
        X = Diagonal(1:3N) # position operator in coordinate representation
        for iφ in eachindex(φₓ)
            XE = tbh.c[:, targetlevels, iφ]' * X * tbh.c[:, targetlevels, iφ] # position operator in energy representation
            pos[:, iφ], d[:, :, iφ] = eigen(Hermitian(XE))
            E[:, iφ] = transpose(tbh.E[targetlevels, iφ]) * abs2.(d[:, :, iφ])
        end
    else
        if periodicity == :space
            X = Diagonal(cis.(2π/(3N*n_t) .* (1:3N*n_t))) # position operator in coordinate representation
        else # if periodicity == :time || periodicity == :both
            X = Diagonal([cis(2π*((n % n_t)/n_t + (n ÷ n_t)/(3N*n_t))) for n in 0:3N*n_t-1]) # position operator in coordinate representation
        end    
        for iφ in eachindex(φₓ)
            XE = tbh.c[:, targetlevels, iφ]' * X * tbh.c[:, targetlevels, iφ] # position operator in energy representation
            _, d[:, :, iφ], pos_complex = schur(XE)
            pos_real = @. mod2pi(angle(pos_complex)) / 2π * 3N*n_t # shift angle from [-π, π) to [0, 2π)
            sp = sortperm(pos_real)                    # sort the eigenvalues
            pos[:, iφ] = pos_real[sp]
            @views Base.permutecols!!(d[:, :, iφ], sp) # sort the eigenvectors in the same way
            E[:, iφ] = transpose(tbh.E[targetlevels, iφ]) * abs2.(d[:, :, iφ])
        end
    end
end

"Calculate Wannier centres for the specified levels of the given Hamiltonian matrix. `n_s` is the number of sites per cell."
function compute_wanniers(H::Array{<:Number, 3}; targetlevels::AbstractVector{<:Integer}, n_s::Integer, isperiodic::Bool)
    N = length(targetlevels)
    pos = Matrix{Float64}(undef, N, size(H, 3))
    if isperiodic
        X = Diagonal(cis.(2π/(n_s*N) * (0:n_s*N-1)))
        for iφ in axes(H, 3)
            c = eigvecs(Hermitian(H[:, :, iφ]))
            XE = c[:, targetlevels]' * X * c[:, targetlevels]
            _, _, pos_complex = schur(XE)
            pos[:, iφ] = (@. mod2pi(angle(pos_complex)) / 2π * n_s*N) |> sort # shift angle from [-π, π) to [0, 2π)
        end
    else
        X = Diagonal(0:n_s*N-1)
        for iφ in axes(H, 3)
            c = eigvecs(Hermitian(H[:, :, iφ]))
            XE = c[:, targetlevels]' * X * c[:, targetlevels]
            pos[:, iφ] = eigvals(Hermitian(XE))
        end
    end
    return pos
end

"""
Construct Wannier functions at each phase number in `whichphases`. All Wannier functions contained in `thb` are constructed.
Return `w`, where `w[:, :, j, i]` = `j`th Wannier function at `whichphases[i]`th phase in the form of a 2D map (first index is temporal, second is spatial).
"""
function make_wannierfunctions(tbh::TBFloquetHamiltonian, whichphases::AbstractVector{<:Integer})
    (;N) = tbh
    targetlevels = tbh.w.targetsubbands # the `targetsubbands` stored in `tbh` are actually level numbers
    n_w = length(targetlevels)
    n_t = 4 # number of temporal sites
    w = Array{ComplexF64, 4}(undef, n_t, 3N, n_w, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for j in 1:n_w
            w[:, :, j, i] = reshape(sum(tbh.w.d[k, j, iφ] * tbh.c[:, targetlevels[k], iφ] for k = 1:n_w), (n_t, 3N))
        end
    end
    return w
end

"Return separated spatial TB Hamiltonian. Most useful when `tbh` describes space-only pumping."
function separate_space(tbh::TBFloquetHamiltonian; N::Integer, periodic::Bool)
    n_t = 4 # total number of temporal sites
    n_s = 3 # number of spatial sites per cell
    H_S = Array{Float64, 3}(undef, n_s*N, n_s*N, size(tbh.H, 3))
    for (iφ, H) in enumerate(eachslice(tbh.H, dims=3))
        J = fill(abs(H[1, n_t+1]), n_s*N-1)
        d = repeat(diag(H)[1:n_t:2n_t+1], N) |> real
        H_S[:, :, iφ] = diagm(0 => d, -1 => J, 1 => J)
        periodic && (H_S[1, end, iφ] = H_S[end, 1, iφ] = J[1])
    end
    return H_S
end

"""
Return separated temporal TB Hamiltonian. Most useful when `tbh` describes time-only pumping.
The method is experimental because it's unclear how to repeat the pattern with phases.
"""
function separate_time(tbh::TBFloquetHamiltonian; N::Integer, periodic::Bool)
    n_t = 2 # number of temporal sites per cell
    H_T = Array{Float64, 3}(undef, n_t*N, n_t*N, size(tbh.H, 3))
    for (iφ, H) in enumerate(eachslice(tbh.H, dims=3))
        J = [repeat([H[1, 2], H[2, 3]], N-1); H[1, 2]]
        d = repeat([H[1, 1], H[2, 2]], N) |> real
        H_T[:, :, iφ] = diagm(0 => d, -1 => conj.(J), 1 => J)
        periodic && (H_T[1, end, iφ] = J[2]; H_T[end, 1, iφ] = J[2]')
    end
    return H_T
end

"Bring the wannier functions contained in `w` to correct order for dispalying animation."
function order_wannierfunctions!(w::Array{Float64, 4}, whichphases::AbstractVector{<:Integer})
    n = size(w, 2) # total number of spatial sites
    n_w = size(w, 3)
    for i in 2:length(whichphases)
        for j in 1:n_w
            prev_maxsite = argmax(w[:, :, j, i-1]) # find which site of the `j`th state contained the maximum at the previous phase
            for k in j+1:n_w
                maxsite = argmax(w[:, :, k, i]) # find which site of the `k`th state contains the maximum at the current phase
                # swap wanniers if we find a state with the same maximum position, or a state whose maximum is shifted according to the pumping protocol
                if prev_maxsite == maxsite ||
                   maxsite - prev_maxsite == CartesianIndex(-1, -1) || maxsite - prev_maxsite == CartesianIndex(1, -1) ||
                   maxsite - prev_maxsite == CartesianIndex(-1, n-1) || maxsite - prev_maxsite == CartesianIndex(1, n-1)
                    temp = w[:, :, j, i]
                    w[:, :, j, i] = w[:, :, k, i]
                    w[:, :, k, i] = temp
                    break
                end
            end
        end
    end
end

export Wanniers,
    UnperturbedHamiltonian,
    diagonalise!,
    make_eigenfunctions,
    compute_wanniers!,
    compute_wanniers,
    make_wannierfunctions,
    SimpleTBHamiltonian,
    kspace_hamiltonian,
    diagonalise_kspace,
    TBHamiltonian,
    FloquetWanniers,
    swap_wanniers!,
    FloquetHamiltonian,
    swap_eigenstates!,
    optimise_wanniers!,
    order_wanniers!,
    TBFloquetHamiltonian,
    separate_space,
    separate_time,
    order_wannierfunctions!
    
end