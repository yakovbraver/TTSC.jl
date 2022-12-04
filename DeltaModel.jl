module DeltaModel

using ProgressMeter: @showprogress
import IntervalRootFinding as iroots
using IntervalArithmetic: (..)
using LinearAlgebra: eigvals, eigen, schur, ⋅, svd, diagm, Diagonal

"A type for storing the Wannier functions."
mutable struct Wanniers
    targetband::Int
    E::Array{Float64, 3} # `E[j, b, i]` = mean energy of `j`th wannier of the `b`th subband (1≤b≤3) at `i`th phase
    pos::Array{Float64, 3} # `pos[j, b, i]` = position eigenvalue of `j`th wannier of the `b`th subband (1≤b≤3) at `i`th phase
    d::Array{ComplexF64, 4} # `d[:, :, b, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details4
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
function diagonalise!(uh::UnperturbedHamiltonian, n_subbands::Integer, bounds::Tuple{<:Real, <:Real})
    (;N, a, U, φₓ) = uh
    uh.E = Array{Float64, 3}(undef, N, n_subbands, length(φₓ)) # the number of different 𝑘's is equal to `N`
    uh.c = Array{ComplexF64, 4}(undef, 6, N, n_subbands, length(φₓ))
    uh.κ = Array{Float64, 4}(undef, 3, N, n_subbands, length(φₓ))

    @showprogress for (iφ, φ) in enumerate(φₓ)
        for ik = 1:N
            ka = 2π*(ik-1)/N

            if ik <= N÷2 + 1
                rts = iroots.roots(ε -> cos_ka(ε; φ, uh) - cos(ka), bounds[1]..bounds[2])
                uh.E[ik, :, iφ] = sort([rts[i].interval.lo for i in eachindex(rts)])
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
                c ./= √(N*X) 
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

"Calculate Wannier vectors for the unperturbed Hamiltonian `uh`."
function compute_wanniers!(uh::UnperturbedHamiltonian, targetband::Integer)
    (;N, a, φₓ, c, E, κ) = uh
    uh.w.targetband = targetband

    X = Matrix{ComplexF64}(undef, N, N) # position operator
    
    k₂ = 2π/(N*a)
    𝐹(x, i, n, j′, j, m, iφ) = begin
        κʲ = κ[i, j, m, iφ]
        κʲ′ = κ[i, j′, m, iφ]
        cis((n-1)*2π*(j-j′)/N - (κʲ′ + κʲ - k₂)x) / 4 * (
            im * (c[2i-1, j′, m, iφ] + im*c[2i, j′, m, iφ])' * (c[2i-1, j, m, iφ] - im*c[2i, j, m, iφ]) / (κʲ′ + κʲ - k₂) +
            (c[2i-1, j, m, iφ] + im*c[2i, j, m, iφ]) * cis(2κʲ*x) * ( 
                (c[2i, j′, m, iφ] - im*c[2i-1, j′, m, iφ]) / (-κʲ′ + κʲ + k₂) +
                (c[2i, j′, m, iφ] + im*c[2i-1, j′, m, iφ]) / ( κʲ′ + κʲ + k₂) * cis(-2κʲ′*x) )' +
            (c[2i-1, j′, m, iφ] - im*c[2i, j′, m, iφ])' * (c[2i, j, m, iφ] + im*c[2i-1, j, m, iφ]) * cis(2κʲ′*x) / (κʲ′ - κʲ + k₂) )
    end
    for iφ in eachindex(φₓ)
        for b in 1:3 # for each of the 3 subbands in the target band
            m = 3(targetband-1) + b # "global" subband number
            for j in 1:N
                for j′ in 1:N
                    X[j′, j] = 0
                    for n = 1:N, i = 1:3
                        X[j′, j] += 𝐹((n-1)a + i*a/3, i, n, j′, j, m, iφ) - 𝐹((n-1)a + (i-1)a/3, i, n, j′, j, m, iφ)
                    end
                end
            end
            # `eigen` does not guarantee orthogonality of eigenvectors in case of degeneracies for `X` unitary, so use `schur`
            # (although a degeneracy of coordinates eigenvalues is unlikely here)
            _, uh.w.d[:, :, b, iφ], pos_complex = schur(X)
            pos_real = @. (angle(pos_complex) + pi) / k₂ # shift angle from [-π, π) to [0, 2π)
            sp = sortperm(pos_real)                         # sort the eigenvalues
            uh.w.pos[:, b, iφ] = pos_real[sp]
            @views Base.permutecols!!(uh.w.d[:, :, b, iφ], sp) # sort the eigenvectors in the same way
            uh.w.E[:, b, iφ] = [abs2.(dˣ) ⋅ E[:, m, iφ] for dˣ in eachcol(uh.w.d[:, :, b, iφ])]
        end
    end
end

"""
Calculate Wannier vectors for the unperturbed Hamiltonians `h` by mixing the states corresponding to 𝑘 = 0 in each of the 3 subbands of the `targetband`.
Return `d, pos, E` as contained in `Wanniers` struct, except that these do not contain a separate dimension for the different 𝑘's. 
"""
function compute_wanniers(uh::UnperturbedHamiltonian, targetband::Integer)
    (;N, a, c, E, κ) = uh

    X = Matrix{ComplexF64}(undef, 3, 3) # position operator
    
    k₂ = 2π/a
    iφ = 1
    ik = 1 # which value of 𝑘 to take in the `c` array; take the first since there is only one
    𝐹(x, i, j′, j) = begin
        κʲ = κ[i, ik, j, iφ]
        κʲ′ = κ[i, ik, j′, iφ]
        cis(-(κʲ′ + κʲ - k₂)x) / 4 * (
            im * (c[2i-1, ik, j′, iφ] + im*c[2i, ik, j′, iφ])' * (c[2i-1, ik, j, iφ] - im*c[2i, ik, j, iφ]) / (κʲ′ + κʲ - k₂) +
            (c[2i-1, ik, j, iφ] + im*c[2i, ik, j, iφ]) * cis(2κʲ*x) * ( 
                (c[2i, ik, j′, iφ] - im*c[2i-1, ik, j′, iφ]) / (-κʲ′ + κʲ + k₂) +
                (c[2i, ik, j′, iφ] + im*c[2i-1, ik, j′, iφ]) / ( κʲ′ + κʲ + k₂) * cis(-2κʲ′*x) )' +
            (c[2i-1, ik, j′, iφ] - im*c[2i, ik, j′, iφ])' * (c[2i, ik, j, iφ] + im*c[2i-1, ik, j, iφ]) * cis(2κʲ′*x) / (κʲ′ - κʲ + k₂) )
    end

    for b in 1:3  # `b` and `b′` run over the 3 subbands
        m = 3(targetband-1) + b # "global" subband number
        for b′ in 1:3
            m′ = 3(targetband-1) + b′
            X[b′, b] = 0
            for i = 1:3
                X[b′, b] += 𝐹(i*a/3, i, m′, m) - 𝐹((i-1)a/3, i, m′, m)
            end
        end
    end
    _, d, pos_complex = schur(X)
    pos_real = @. (angle(pos_complex) + pi) / k₂ # shift angle from [-π, π) to [0, 2π)
    sp = sortperm(pos_real)          # sort the eigenvalues
    pos = pos_real[sp]
    @views Base.permutecols!!(d, sp) # sort the eigenvectors in the same way
    E = [sum(abs2(dˣ[b]) * uh.E[ik, 3(targetband-1) + b, iφ] for b in 1:3) for dˣ in eachcol(d)]
    return d, pos, E
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
            for j in 1:uh.N
                w[:, j, b, i] = sum(uh.w.d[p, j, b, iφ] * ψ[:, p, b, i] for p = 1:uh.N)
            end
        end
    end
    return x, ψ, w
end

"""
A type representing the tight-binding Hamiltonian
    hₜ = ∑ⱼ 𝐽₁𝑏⁺ⱼ𝑎ⱼ + 𝐽₂𝑐⁺ⱼ𝑏ⱼ + 𝐽₃𝑎⁺ⱼ₊₁𝑐ⱼ + h.c.
         + 𝑈 ∑ⱼ 𝑎⁺ⱼ𝑎ⱼcos(𝜑ₓ) + 𝑏⁺ⱼ𝑏ⱼcos(𝜑ₓ + 2π/3) + 𝑐⁺ⱼ𝑐ⱼcos(𝜑ₓ + 4π/3)
"""
mutable struct TBHamiltonian
    N::Int # number of lattice cells
    a::Float64
    U::Float64
    J::Vector{ComplexF64}
    isperiodic::Bool
    φₓ::Vector{Float64}
    E::Matrix{Float64}      # `E[i, j]` = `i`th eigenvalue at `j`th phase, `i` ∈ [1, `3N`], `j` ∈ [1, `length(φₓ)`]
    c::Array{ComplexF64, 3} # `c[:, i, j]` = `i`th eigenvector at `j`th phase
    w::Wanniers 
end

"Construct a `TBHamiltonian` object."
function TBHamiltonian(n_cells::Integer; a::Real, U::Real, J::Vector{<:Number}, isperiodic::Bool, φₓ::AbstractVector{<:Real})
    E = Matrix{Float64}(undef, 3n_cells, length(φₓ))
    c = Array{ComplexF64, 3}(undef, 3n_cells, 3n_cells, length(φₓ))
    w = Wanniers(0, Array{Float64,3}(undef, n_cells, 3, length(φₓ)), Array{Float64,3}(undef, n_cells, 3, length(φₓ)),
                 Array{ComplexF64,4}(undef, n_cells, n_cells, 3, length(φₓ)))
    TBHamiltonian(Int(n_cells), Float64(a), Float64(U), ComplexF64.(J), isperiodic, collect(Float64, φₓ), E, c, w)
end

"Diagonalise the TB Hamiltonian `h` at each phase."
function diagonalise!(h::TBHamiltonian)
    (;N, U, J, isperiodic, φₓ) = h
    for (i, φ) in enumerate(φₓ)
        diag = repeat([U*cos(φ), U*cos(φ+2π/3), U*cos(φ+4π/3)], N)
        J_diag = [repeat(J, N-1); J[1:2]]
        H = diagm(0 => diag, -1 => J_diag, 1 => conj.(J_diag))
        if isperiodic
            H[1, end] = J[3]
            H[end, 1] = J[3]'
        end
        h.E[:, i], h.c[:, :, i] = eigen(H)
    end
end

"Calculate Wannier vectors for each of the three subbands for the TB Hamiltonian `h`."
function compute_wanniers!(h::TBHamiltonian)
    (;N, a, φₓ) = h
    for b in 1:3
        levels = N*(b-1)+1:N*b
        if h.isperiodic
            X = Diagonal([cis(2π/(N*a) * n*a/3) for n in 0:3N-1]) # position operator in coordinate representation
            for iφ in eachindex(φₓ)
                XE = h.c[:, levels, iφ]' * X * h.c[:, levels, iφ] # position operator in energy representation
                _, h.w.d[:, :, b, iφ], pos_complex = schur(XE)
                pos_real = @. mod2pi(angle(pos_complex)) / 2π * N*a # shift angle from [-π, π) to [0, 2π)
                sp = sortperm(pos_real)                        # sort the eigenvalues
                h.w.pos[:, b, iφ] = pos_real[sp]
                @views Base.permutecols!!(h.w.d[:, :, b, iφ], sp) # sort the eigenvectors in the same way
                h.w.E[:, b, iφ] = [abs2.(dˣ) ⋅ h.E[levels, iφ] for dˣ in eachcol(h.w.d[:, :, b, iφ])]
            end
        else
            X = Diagonal([n*a/3 for n in 0:3N-1]) # position operator in coordinate representation
            for iφ in eachindex(φₓ)
                XE = h.c[:, levels, iφ]' * X * h.c[:, levels, iφ] # position operator in energy representation
                h.w.pos[b, :, iφ], h.w.d[:, :, b, iφ] = eigen(XE)
                h.w.E[b, :, iφ] = [abs2.(dˣ) ⋅ h.E[levels, iφ] for dˣ in eachcol(h.w.d[:, :, b, iφ])]
            end
        end
    end
end

"""
Construct Wannier functions at each phase number in `whichphases`. All Wannier functions contained in `uh` are constructed.
Return `w`, where `w[:, j, b i]` = `j`th Wannier function of `b`th subband at `i`th phase.
"""
function make_wannierfunctions(h::TBHamiltonian, whichphases::AbstractVector{<:Integer})
    (;N) = h
    w = Array{ComplexF64, 4}(undef, size(h.c, 1), N, 3, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for b in 1:3
            for j in 1:N
                w[:, j, b, i] = sum(h.w.d[k, j, b, iφ] * h.c[:, N*(b-1)+k, iφ] for k = 1:N)
            end
        end
    end
    return w
end

"Return the 𝑘-space Hamiltonian matrix for `h` at the given phase `φ` and at 𝑘𝑎 = `ka`."
function kspace_hamiltonian(h::TBHamiltonian, φ::Real, ka::Real)
    (;U, J) = h
    [U*cos(φ)        J[1]'         J[3]cis(-ka)
     J[1]            U*cos(φ+2π/3) J[2]'
     (J[3]cis(-ka))' J[2]          U*cos(φ+4π/3)]
end

"""
Diagonalise the TB Hamiltonian `h` in 𝑘-space at each phase for the values of 𝑘𝑎 in `ka`.
Return the matrix of eigenenergies `E`, where `E[:, i]` is the energy at `i`th phase.
In `E`, rows 1:3 corresopnd to `ka[1]`, rows 4:6 correspond to `ka[2]`, and so on.
"""
function diagonalise_kspace(h::TBHamiltonian, ka::AbstractVector{<:Real})
    E = Matrix{Float64}(undef, 3length(ka), length(h.φₓ))
    for (iφ, φ) in enumerate(h.φₓ)
        for ik in eachindex(ka)
            E[3(ik-1)+1:3ik, iφ] .= eigvals(kspace_hamiltonian(h, φ, ka[ik]))
        end
    end
    return E
end

"A type for storing the Wannier functions."
mutable struct FloquetWanniers
    targetsubbands::Vector{Int}
    E::Array{Float64, 2} # `E[j, i]` = mean energy of `j`th wannier at `i`th phase
    pos::Array{Float64, 2} # `pos[j, i]` = position eigenvalue of `j`th wannier at `i`th phase
    d::Array{ComplexF64, 3} # `d[:, :, i]` = position eigenvectors at `i`th phase; see methods of `compute_wanniers!` for details4
end

"Default-construct an empty `FloquetWanniers` object."
FloquetWanniers() = FloquetWanniers(Int[], Float64[;;], Float64[;;], ComplexF64[;;;])

"""
A type representing the Floquet Hamiltonian
    ℋ = ℎ - i∂ₜ + λₛcos²(2𝑥)cos(2𝜔𝑡) + λₗcos²(2𝑥)cos(𝜔𝑡 + 𝜑ₜ),
where ℎ is the unperturbed Hamiltonian represented by [`UnperturbedHamiltonian`](@ref), and 𝜑ₜ = 2𝜑ₓ.
"""
mutable struct FloquetHamiltonian
    uh::UnperturbedHamiltonian
    s::Int
    λₛ::Float64
    λₗ::Float64
    ω::Float64
    pumptype::Symbol
    E::Array{Float64, 3} # `E[i, ik, j]` = `i`th eigenvalue (Floquet quasienergy) at `j`th phase, `i = 1:maxlevel`
    b::Array{ComplexF64, 4} # `b[:, i, ik, j]` = `i`th eigenvector at `j`th phase, `i = 1:maxlevel; j = 1:2maxlevel+1`
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

"Diagonalise the Floquet Hamiltonian `fh` at each phase."
function diagonalise!(fh::FloquetHamiltonian)
    (;N, a, φₓ, E, c, κ) = fh.uh
    (;s, ω, λₛ, λₗ, pumptype, ν) = fh

    n_levels = size(fh.E, 1)

    H = zeros(ComplexF64, n_levels, n_levels) # ℋ matrix

    𝐹(x, i, ik, j′, j, iφ, k₂) = begin
        κʲ = κ[i, ik, j, iφ]
        κʲ′ = κ[i, ik, j′, iφ]
        cis(-(κʲ′ + κʲ - k₂)x) / 4 * (
            im * (c[2i-1, ik, j′, iφ] + im*c[2i, ik, j′, iφ])' * (c[2i-1, ik, j, iφ] - im*c[2i, ik, j, iφ]) / (κʲ′ + κʲ - k₂) +
            (c[2i-1, ik, j, iφ] + im*c[2i, ik, j, iφ]) * cis(2κʲ*x) * ( 
                (c[2i, ik, j′, iφ] - im*c[2i-1, ik, j′, iφ]) / (-κʲ′ + κʲ + k₂) +
                (c[2i, ik, j′, iφ] + im*c[2i-1, ik, j′, iφ]) / ( κʲ′ + κʲ + k₂) * cis(-2κʲ′*x) )' +
            (c[2i-1, ik, j′, iφ] - im*c[2i, ik, j′, iφ])' * (c[2i, ik, j, iφ] + im*c[2i-1, ik, j, iφ]) * cis(2κʲ′*x) / (κʲ′ - κʲ + k₂) )
    end

    for ik in 1:N
        for (iφ, φ) in enumerate(φₓ)
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
                        for i = 1:3, k₂ in (-4, 4)
                            ∫cos += 𝐹(i*a/3, i, ik, m′, m, iφ, k₂) - 𝐹((i-1)a/3, i, ik, m′, m, iφ, k₂)
                        end
                        # if pumping is space-time, then also multiply by cis(-𝜑ₜ). `φ` runs over 𝜑ₓ, and we assume the pumping protocol 𝜑ₜ = 𝜑ₓ
                        H[m′, m] = (pumptype == :space ? λₗ/8 * ∫cos : λₗ/8 * ∫cos * cis(-φ))
                    elseif pumptype == :time 
                        H[m′, m] *= cis(-(φₓ[iφ]-φₓ[iφ-1]))
                    end
                    H[m, m′] = H[m′, m]'
                end
                
                # place the elements of the short lattice
                for g in 1:3
                    m′ = 3(2s + ν[m] - 1) + g
                    m′ > n_levels && break
                    if pumptype != :time || iφ == 1 # if pumping is time-only, this must be calculated only once, at `iφ` = 1
                        ∫cos = ComplexF64(0)
                        for i = 1:3, k₂ in (-4, 4)
                            ∫cos += 𝐹(i*a/3, i, ik, m′, m, iφ, k₂) - 𝐹((i-1)a/3, i, ik, m′, m, iφ, k₂)
                        end
                        H[m′, m] = λₛ/8 * ∫cos
                    end
                    H[m, m′] = H[m′, m]'
                end
            end
            fh.E[:, ik, iφ], fh.b[:, :, ik, iφ] = eigen(H)
        end
    end
end

"""
Permute Floquet quasienergy levels contained in `fh.E` so that they are stored in the same order as the eigenenergies of ℎ stored in `fh.uh.E`.
Repeat this for every phase.
To perfrorm the sorting, first calculate `fh.uh.E - fh.ν[m]`, which is the diagonal of ℋ. If there is no perturbation, then these
are the Floquet quasienergies. Then, sort them in ascending order (as if we diagonalised the Hamiltonian) and find the permutation
that would undo this sorting. This permutation is applied to a copy of `fh.E`.
The procedure yields fully correct results only if `fh.E` has been calculated at zero perturbation. The perturbation may additionally change
the order of levels, and there is no simple way of disentangling the order. The permutation is still useful in that case, but the results 
should not be taken too literally.
"""
function order_floquet_levels(fh::FloquetHamiltonian)
    E = similar(fh.E)
    for iφ in axes(fh.E, 3)
        for ik in axes(fh.E, 2)
            E_diag = [fh.uh.E[ik, m, iφ] - fh.ν[m] * fh.ω/fh.s for m in axes(fh.uh.E, 2)] # Floquet energies at zero perturbation
            invsort = sortperm(sortperm(E_diag))  # inverse permutation, such that `sort(E_diag)[invsort] == E_diag`
            E[:, ik, iφ] .= fh.E[invsort, ik, iφ]
        end
    end
    return E
end

"""
Construct Floquet modes at coordinates in `x` and time moments in `Ωt` for each spatial subband number in `whichsubband` at each phase number in `whichphases`.
Return `x, u`, where `x` are the abscissas and `u[ix, it, ik, j, i]` = wavefunction corresponding to the `ik`th value of 𝑘 and `j`th spatial subband at `i`th phase at `ix`th coordinate at `it`th time moment.
"""
function make_eigenfunctions(fh::FloquetHamiltonian, n_x::Integer, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer},
                             whichsubbands::AbstractVector{<:Integer})
    (;N, a) = fh.uh
    x = range(0, a*N, 3N*n_x+1)
    u = Array{ComplexF64, 5}(undef, 3N*n_x+1, length(Ωt), N, length(whichsubbands), length(whichphases))
    n_levels = size(fh.E, 1) # number of Floquet levels; equivalently, number of levels of ℎ used for constructing ℋ
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
"""
function compute_wanniers!(fh::FloquetHamiltonian; targetsubbands::AbstractVector{<:Integer})
    (;N, a, φₓ, c, κ) = fh.uh

    fh.w.targetsubbands = targetsubbands # save this because it's needed in `make_wannierfunctions`

    n_w = length(targetsubbands) * N
    E = Matrix{Float64}(undef, n_w, length(φₓ))
    pos = Matrix{Float64}(undef, n_w, length(φₓ))
    d = Array{ComplexF64, 3}(undef, n_w, n_w, length(φₓ))
    fh.w = FloquetWanniers(targetsubbands, E, pos, d)

    X = Matrix{ComplexF64}(undef, n_w, n_w) # position operator
    
    n_levels = size(fh.E, 1)
    n_subbands = length(targetsubbands)

    expik = Array{ComplexF64, 4}(undef, n_levels, n_levels, N, N)

    k₂ = 2π/(N*a)
    𝐹(x, i, n, j′, j, m′, m, iφ) = begin
        κʲ = κ[i, j, m, iφ]
        κʲ′ = κ[i, j′, m′, iφ]
        cis((n-1)*2π*(j-j′)/N - (κʲ′ + κʲ - k₂)x) / 4 * (
            im * (c[2i-1, j′, m′, iφ] + im*c[2i, j′, m′, iφ])' * (c[2i-1, j, m, iφ] - im*c[2i, j, m, iφ]) / (κʲ′ + κʲ - k₂) +
            (c[2i-1, j, m, iφ] + im*c[2i, j, m, iφ]) * cis(2κʲ*x) * ( 
                (c[2i, j′, m′, iφ] - im*c[2i-1, j′, m′, iφ]) / (-κʲ′ + κʲ + k₂) +
                (c[2i, j′, m′, iφ] + im*c[2i-1, j′, m′, iφ]) / ( κʲ′ + κʲ + k₂) * cis(-2κʲ′*x) )' +
            (c[2i-1, j′, m′, iφ] - im*c[2i, j′, m′, iφ])' * (c[2i, j, m, iφ] + im*c[2i-1, j, m, iφ]) * cis(2κʲ′*x) / (κʲ′ - κʲ + k₂) )
    end

    for iφ in eachindex(φₓ)
        # if pumping is time-only, then `expik` must be calculated only at the first iteration, thereby using `c`'s at 𝜑ₓ = 0
        if fh.pumptype != :time || iφ == 1
            for ik in 1:N,  ik′ in 1:N,  m in 1:n_levels,  m′ in 1:n_levels
                expik[m′, m, ik′, ik] = 0
                for n = 1:N, i = 1:3
                    expik[m′, m, ik′, ik] += 𝐹((n-1)a + i*a/3, i, n, ik′, ik, m′, m, iφ) - 𝐹((n-1)a + (i-1)a/3, i, n, ik′, ik, m′, m, iφ)
                end
            end
        end

        t = (fh.pumptype == :space ? π/5 : π/5 - iφ/length(φₓ)*π/2) # time moment at which to diagonalise the coordinate operator

        for ik in 1:N,  ik′ in 1:N
            for (in, n) in enumerate(targetsubbands)
                for (in′, n′) in enumerate(targetsubbands)
                    X[in′+(ik′-1)*n_subbands, in+(ik-1)*n_subbands] = sum(fh.b[m, n, ik, iφ] * sum(fh.b[m′, n′, ik′, iφ]' * expik[m′, m, ik′, ik] *
                                                                          cis((fh.ν[m′] - fh.ν[m]) * t) for m′ in 1:n_levels) for m in 1:n_levels)
                end
            end
        end
        _, d[:, :, iφ], pos_complex = schur(X)
        pos_real = @. (angle(pos_complex) + pi) / k₂ # shift angle from [-π, π) to [0, 2π)
        sp = sortperm(pos_real)         # sort the eigenvalues
        pos[:, iφ] = pos_real[sp]
        @views Base.permutecols!!(d[:, :, iφ], sp) # sort the eigenvectors in the same way
        E[:, iφ] = [sum(abs2(dˣ[m]) * fh.E[targetsubbands[(m-1)%n_subbands+1], (m-1)÷n_subbands+1, iφ] for m in eachindex(dˣ)) for dˣ in eachcol(d[:, :, iφ])]
    end
end

"""
Construct Wannier functions at coordinates in `x` at each phase number in `whichphases`. All Wannier functions contained in `fh` are constructed.
In the process, energy eigenfunctions are also constructed.
Return `u, w`, where `w[ix, it, j, i]` = `j`th Wannier function at `i`th phase at `ix`th coordinate at `it`th time moment,
and `u` is an array of Floquet modes in the same format.
"""
function make_wannierfunctions(fh::FloquetHamiltonian, n_x::Integer, Ωt::AbstractVector{<:Real}, whichphases::AbstractVector{<:Integer})
    (;N) = fh.uh
    n_subbands = length(fh.w.targetsubbands)
    n_w = n_subbands * N
    x, u = make_eigenfunctions(fh, n_x, Ωt, whichphases, fh.w.targetsubbands) # `u[ix, it, ik, j, i]`
    w = Array{ComplexF64, 4}(undef, length(x), length(Ωt), n_w, length(whichphases))
    for (i, iφ) in enumerate(whichphases)
        for j in 1:n_w
            w[:, :, j, i] = sum(fh.w.d[m, j, iφ] * u[:, :, (m-1)÷n_subbands+1, (m-1)%n_subbands+1, i] for m = 1:n_w)
        end
    end
    return x, u, w
end

end