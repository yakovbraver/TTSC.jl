using Plots, UnPack, LaTeXStrings, ProgressBars
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

"""
Hamiltonian (S17) defined here as
    𝐻 = (ℎ₀ - Ω𝐼₀) + (𝐼 - 𝐼₀)²/2𝑀 + λL⋅cos(𝑠ϑ + χL) + λS⋅cos(2𝑠ϑ + χS)
"""
Base.@kwdef mutable struct Heff
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
function getE(H::Heff, ϑ::AbstractVector, I::AbstractVector)
    @unpack h0, Ω, I0, M, λL, χL, λS, χS, s = H
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = h0 - Ω*I0 + (I[i]-I0)^2/2M + λL*cos(s*ϑ[t] + χL) + λS*cos(2s*ϑ[t] + χS)
    end
    E
end

function plot_isoenergies(H::Heff)
    ϑ = range(0, 2π, length=50)
    I = range(H.I0-10, H.I0+10, length=50)
    E = getE(H, ϑ, I)
    # calculate fixed energy limits so that the z-scale is the same regardless of the phases
    Elims = H.h0 - H.Ω*H.I0 .+ (-H.λL - H.λS, H.λL + H.λS + (I[end]-H.I0)^2/2H.M)
    fig = contourf(ϑ, I, E', xlabel=L"\theta"*", rad", ylabel=L"I", cbartitle="Energy", levels=range(Elims[1], Elims[2], length=20), color=:viridis)
    title!(L"M = %$(round(H.M, sigdigits=2)),"*
           L"\lambda_L = %$(round(H.λL, sigdigits=2)), \chi_L = %$(round(H.χL, sigdigits=2)),"*
           L"\lambda_S = %$(round(H.λS, sigdigits=2)), \chi_S = %$(round(H.χS, sigdigits=2))")
    fig
end

# e = plot_isoenergies(Heff(h0=2.0, Ω=2.0, I0=1.0, M=0.5, λL=0.01, χL=0.2, λS=0.02, χS=0.0, s=2))
ϕ = range(0, 2π, length=50)
H = Heff(h0=0.0, Ω=0.0, I0=0.0, M=0.5, λL=10, χL=0, λS=10, χS=0.0, s=2)
@gif for i in ProgressBar(ϕ)
    H.χL = i
    plot_isoenergies(H)
end

using BandedMatrices, KrylovKit

function plot_bands(; M, λL, λS)
    N_eigvals = 10 # number of eigenvalues to compute for each 𝑘
    N_j = N_eigvals # number of indices 𝑗 to take for the Hamiltonian (its size will be (2N_j+1)×(2N_j+1))
    
    # Hamiltonian matrix
    H = BandedMatrix{ComplexF64}(undef, (2N_j + 1, 2N_j + 1), (2, 2))
    H[band(-2)] .= repeat([λS], 2N_j-1)
    H[band(2)] .= repeat([λS], 2N_j-1)
    
    N_ϕ = 20 # number of adiabatic phases to sweep over
    # `bands` stores the dalculated eigenvalues: columns numerate the adiabatic phase, while rows numerate eigenvalues.
    # Rows `1:N_eigvals` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
    # Rows `N_eigvals:end` store the eigenvalues corresponding to the boundary of BZ, in our case 𝑘 = 1.
    bands = Matrix{Float64}(undef, 2N_eigvals, N_ϕ)

    phases = range(0, 2π, length=N_ϕ)
    for k in [0, 1]
        H[band(0)] .= [(2j + k)^2 / M for j = -N_j:N_j]
        a = k*N_eigvals + 1 # controls where to place the eigenvalues depedning on `k`; see description of `bands`
        b = a+N_eigvals - 1 # controls where to place the eigenvalues depedning on `k`; see description of `bands`
        for (i, ϕ) in enumerate(phases)
            H[band(-1)] .= repeat([λL*cis(-ϕ)], 2N_j)
            H[band(1)] .= repeat([λL*cis(ϕ)], 2N_j)
            vals, _, _ = eigsolve(H, N_eigvals, :SR)
            bands[a:b, i] .= vals[1:N_eigvals] ./ 4π # we do not include the factor 4π in the diagonalisation problem and restore it here
        end
    end
    fig = plot()
    for i in 1:N_eigvals
        plot!(fig, phases, bands[i, :], fillrange=bands[N_eigvals+i, :], fillalpha=0.35, label="band $i")
    end
    xlabel!(L"\varphi_t"*", rad")
    ylabel!("Energy")
    fig
end

M = 0.5; λL = 10; λS = 10;
f = plot_bands(; M, λL, λS)

# N = 20
# d = [repeat([2, 1], N÷2-1); 2]
# ssh = BandedMatrix(-1 => d, 1 => d)
# vals, vecs, info = eigsolve(ssh, N, :SR, tol=1e-20)
# scatter(repeat([2], N), vals)