using Plots, LaTeXStrings
using SpecialFunctions: gamma
using Combinatorics: factorial, binomial
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

### Calculate the spatial Hamiltonian as a function of the action, and compute the associated derivatives

include("SpacetimeHamiltonian.jl")

function 𝐻₀(p, x, params)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2
end

function 𝐻(p, x, params, t)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2 +
    params[4]*sin(2x)^2*cos(2params[6]*t) + 
    params[5]*cos(2x)^2*cos(params[6]*t)
end

function 𝑄ₛ(p::Real, x::Real)
    sin(2x)^2
end

function 𝑄ₗ(p::Real, x::Real)
    cos(2x)^2
end

g = 5000; l = 2;
gₗ = 2g*factorial(l) / √π / gamma(l + 0.5)
Vₗ = 2000

λₛ = 400; λₗ = 50; ω = 540;
s = 2
params = [gₗ, l, Vₗ, λₛ, λₗ, ω]
plot!(range(0, 2π, length=200), x -> 𝐻₀(0, x, params))
H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, (0.8, 1), (1.2, 1.8), 0.05)

function plot_actions(H::SpacetimeHamiltonian)
    figs = [plot() for _ in 1:4];
    x = range(0, π, length=200);
    figs[1] = plot(x, H.𝑈, xlabel=L"x", ylabel=L"U(x)=\tilde{g}_\ell\cos^{2\ell}(2x)+V_L\cos^{2}(x)", legend=false);
    title!(L"\ell = %$l, g = %$g, V_L = %$Vₗ");
    I = Dierckx.get_knots(H.𝐸)
    figs[2] = plot(I, H.𝐸(I), xlabel=L"I", ylabel=L"E", legend=false);
    figs[3] = plot(I, H.𝐸′, xlabel=L"I", ylabel=L"dE/dI", legend=false);
    figs[4] = plot(I, H.𝐸″, xlabel=L"I", ylabel=L"d^2E/dI^2", legend=false, ylims=(-20, 20));
    lay = @layout [a{0.5w} grid(3,1)]
    plot(figs..., layout=lay)
    # savefig("H0.pdf")
end

plot_actions(H)
savefig("H0-parameters.pdf")

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation

function plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, Iₛ, s, levels::Union{Nothing, Vector{<:AbstractFloat}}=nothing)
    ϑ = range(0, 2π, length=50)
    I_max = last(Dierckx.get_knots(H.𝐸))
    I = range(0, I_max, length=50)
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    h₀ = H.𝐸(Iₛ) - ω/s*Iₛ
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] + χₛ) + λₗ*Aₗ*cos(s*ϑ[t] + χₗ)
    end
    levels === nothing ? contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (S17)", color=:viridis) :
                         contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (S17)", color=:viridis; levels)
    hline!([Iₛ], label=L"I_s = %$(round(Iₛ, sigdigits=4))", c=:white)
    title!(L"\omega = %$ω, s = %$s"*"\n"*
    L"\lambda_L = %$(round(λₗ, sigdigits=2)), A_L = %$(round(Aₗ, sigdigits=2)), \chi_L = %$(round(χₗ, sigdigits=2)),"*"\n"*
    L"\lambda_S = %$(round(λₛ, sigdigits=2)), A_S = %$(round(Aₛ, sigdigits=2)), \chi_S = %$(round(χₛ, sigdigits=2))")
end

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [-2s, -s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
ϕₜ = 0.0
eQ = cis(ϕₜ)*coeffs[2]
Aₗ = abs(eQ); χₗ = angle(eQ)

plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, Iₛ, s)
levels = [range(-1000, 500, length=10); range(501, 750, length=10)]
plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, Iₛ, s, levels)
savefig("secular-isoenergies.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

fig = plot();
for i in 1:27
    I, Θ = compute_IΘ(H, i, n_T=200, ϑ₀=0.0)
    scatter!(mod2pi.(Θ.+pi/2), I, xlabel=L"\theta", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
end
for i in [22; 22.5; 23:26]
    I, Θ = compute_IΘ(H, i, n_T=200, ϑ₀=0.75)
    scatter!(mod2pi.(Θ.+pi/2), I, xlabel=L"\theta", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
end
ylims!((0, last(Dierckx.get_knots(H.𝐸))))
title!(L"\ell = %$l, g = %$g, V_L = %$Vₗ, \lambda_S = %$λₛ, \lambda_L = %$λₗ, \omega = %$ω")
display(fig)
savefig("exact-isoenergies.pdf")

### Calculate bands

import BandedMatrices as BM
using SparseArrays: sparse
using KrylovKit: eigsolve

"""
Calculate `n_bands` energy bands of Hamiltonian (S32) sweeping over the adiabatic `phases` (φₜ in (S32)).
In the returned matrix of bands, columns numerate the adiabatic phase, while rows numerate eigenvalues.
Rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
Rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case λₗAₗcos(sϑ+φₜ) leads to 𝑘 = s/2.
"""
function compute_bands(; n_bands::Integer, phases::AbstractVector, s::Integer, M::Real, λₗAₗ::Real, λₛAₛ::Real)
    n_j = 4n_bands # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
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
            vals, _, _ = eigsolve(H, n_bands, :SR)
            bands[a:b, i] .= vals[1:n_bands]
        end
    end
    return bands / 4π # we do not include the factor 4π in the diagonalisation problem and restore it here
end

phases = range(0, 2π, length=40) # values of the adiabatic phase in (S32)
n_bands = 4
bands = compute_bands(; n_bands, phases, s, M=abs(M), λₗAₗ=λₗ*Aₗ, λₛAₛ=λₛ*Aₛ)

fig1 = plot();
for i in 1:n_bands
    plot!(phases, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i");
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Energy")
savefig("bands.pdf")

### Extract tight-binding parameters

gap = bands[2, 1] - bands[1, 1]
w = gap/2 - bands[2, 1]

function tb_parameters(E_0_0, E_0_pi, E_k_pi, k)
    @show J₀ = E_0_0 / 2
    @show Δ = √(E_0_pi^2 - 4J₀^2)
    ε = (E_k_pi^2 - Δ^2 - 2J₀^2 * (1+cos(k))) / (2J₀^2 * (1-cos(k))) |> sqrt
    return J₀, Δ, ε
end

J₀, Δ, ε = tb_parameters(0.235+w, 1.317+w, 1.053+w, 1)
J₀ = 0.23412895701313893
Δ = 1.4778478020259276
ε = 0.5
E0 = @. sqrt(Δ^2*sin(phases)^2 + 4J₀^2)
plot!(phases, E0 .- w)
k = 1
Ek = @. Δ^2*sin(phases)^2 + 2J₀^2 * (1+cos(k) + ε^2*sin(phases)^2 * (1-cos(k))) |> sqrt
plot!(phases, Ek .- w)

### (S23)

"""
Calculate `n_bands` energy bands of Hamiltonian (S20) sweeping over the adiabatic `phases` φₓ and φₜ.
Return a tuple of a matrix `ϵₖ` of `2n_bands` bands of ℎₖ and a matrix `Eₖ` of `n_bands` bands of 𝐻ₖ.
In the returned matrices, columns numerate the adiabatic phases, while rows numerate eigenvalues.
In `Eₖ`, rows `1:n_bands` store the eigenvalues corresponding to the centre of BZ, 𝑘 = 0.
In `Eₖ`, rows `n_bands:end` store the eigenvalues corresponding to the boundary of BZ, in our case Vₗcos²(x+φₓ) leads to 𝑘 = 2/2 = 1.
The dimension of the constructed 𝐻ₖ matrix will be `2n_bands`, hence that many eigenvalues of ℎₖ will be required. This in turn
required constructing ℎₖ of dimension `4n_bands`.
"""
function compute_bands_exact(; n_bands::Integer, phases::AbstractVector, s::Integer, l::Real, gₗ::Real, Vₗ::Real, λₗ::Real, λₛ::Real, ω::Real)
    n_j = 2n_bands  # number of indices 𝑗 to use for constructing the Hamiltonian (its size will be (2n_j+1)×(2n_j+1))
    
    hₖ = BM.BandedMatrix(BM.Zeros{ComplexF64}(2n_j+1, 2n_j+1), (2l, 2l))
    # fill the off-diagonals with binomial numbers; the diagonal is treated in the `k` loop
    for n in 1:l
        hₖ[BM.band(2n)] .= hₖ[BM.band(-2n)] .= gₗ / 4^l * binomial(2l, l-n)
    end
    
    ϵₖ = Matrix{Float64}(undef, 2n_j, length(phases)) # eigenvalues of ℎₖ (eigenenergies of the unperturbed Hamiltonian)
    cₖ = [Vector{ComplexF64}(undef, 2n_j+1) for _ in 1:n_j]  # eigenvectors of ℎₖ
    
    Eₖ = Matrix{Float64}(undef, 2n_bands, length(phases)) # eigenvalues of 𝐻ₖ (Floquet quasi-energies) that will be saved; size is twice `n_bands` for the two values of `k``
    Hₖ_dim = 2n_bands # dimension of the constructed 𝐻ₖ matrix (twice larger than the number of requested quasi-energies)
    n_Hₖ_nonzeros = 9Hₖ_dim - 24s # number of non-zero elements in 𝐻ₖ
    Hₖ_rows = Vector{Int}(undef, n_Hₖ_nonzeros)
    Hₖ_cols = Vector{Int}(undef, n_Hₖ_nonzeros)
    Hₖ_vals = Vector{ComplexF64}(undef, n_Hₖ_nonzeros)
    for k in [0, 1] # iterate over the centre of BZ and then the boundary
        hₖ[BM.band(0)] .= [(2j + k)^2 + Vₗ/2 + gₗ / 4^l * binomial(2l, l) for j = -n_j:n_j]
        # `a` and `b` control where to place the eigenvalues of 𝐻ₖ and ℎₖ depedning on `k`; see function docstring
        a_Hₖ = (k > 0)*n_bands + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_Hₖ = a_Hₖ+n_bands - 1
        a_hₖ = (k > 0)*n_j + 1 # `(k > 0)` is zero for BZ centre (when `k == 0`) and unity otherwise
        b_hₖ = a_hₖ+n_j - 1
        for (z, ϕ) in enumerate(phases)
            hₖ[BM.band(-1)] .= Vₗ/4 * cis(2ϕ)
            hₖ[BM.band(1)]  .= Vₗ/4 * cis(-2ϕ)
            vals, vecs, info = eigsolve(hₖ, n_j, :SR; tol=1.0, krylovdim=n_j+10)
            if info.converged < n_j
                @warn "Only $(info.converged) eigenvalues out of $(n_j) converged when diagonalising ℎₖ. Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
            end
            ϵₖ[a_hₖ:b_hₖ, z] = vals[1:n_j]
            cₖ .= vecs[1:n_j]
            # println(info)

            # Construct 𝐻ₖ
            p = 1 # a counter for placing elements to the vectors Hₖ_*
            for m in 1:Hₖ_dim
                # place the diagonal element (S25)
                Hₖ_rows[p] = Hₖ_cols[p] = m
                Hₖ_vals[p] = ϵₖ[m, z] - ceil(m/2)*ω/s
                p += 1

                # place the elements of the long lattice (S26)
                for i in 1:2
                    m′ = 2s + 2(ceil(Int, m/2)-1) + i
                    m′ > Hₖ_dim && break
                    Hₖ_rows[p] = m′
                    Hₖ_cols[p] = m
                    # the index should run as `j = -n_j+2:n_j-2`, but we don't have negative indexes in the vector, so 
                    j_sum = sum( (cₖ[m′][j+2]/4 + cₖ[m′][j-2]/4 + cₖ[m′][j]/2)' * cₖ[m][j] for j = 3:2n_j-1 ) + 
                            (cₖ[m′][3]/4 + cₖ[m′][1]/2)' * cₖ[m][1] +                # iteration j = 1
                            (cₖ[m′][2n_j-1]/4 + cₖ[m′][2n_j+1]/2)' * cₖ[m][2n_j+1]   # iteration j = 2n_j+1
                    Hₖ_vals[p] = λₗ * cis(-ϕ)/2 * j_sum
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
                    j_sum = sum( (-cₖ[m′][j+2]/4 - cₖ[m′][j-2]/4 + cₖ[m′][j]/2)' * cₖ[m][j] for j = 3:2n_j-1 ) + 
                            (-cₖ[m′][3]/4 + cₖ[m′][1]/2)' * cₖ[m][1] +                # iteration j = 1
                            (-cₖ[m′][2n_j-1]/4 + cₖ[m′][2n_j+1]/2)' * cₖ[m][2n_j+1]   # iteration j = 2n_j+1
                    Hₖ_vals[p] = λₛ/2 * j_sum
                    p += 1
                    # place the conjugate element
                    Hₖ_rows[p] = m
                    Hₖ_cols[p] = m′
                    Hₖ_vals[p] = Hₖ_vals[p-1]'
                    p += 1
                end
            end
            Hₖ = sparse(Hₖ_rows, Hₖ_cols, Hₖ_vals)
            vals, vecs, info = eigsolve(Hₖ, n_bands, :SR; tol=1.0, krylovdim=n_bands+10)
            if info.converged < n_bands
                @warn "Only $(info.converged) eigenvalues out of $(n_bands) converged when diagonalising 𝐻ₖ. Results may be inaccurate." unconverged_norms=info.normres[info.converged+1:end]
            end
            Eₖ[a_Hₖ:b_Hₖ, z] .= vals[1:n_bands]
        end
    end
    # return ϵₖ
    return ϵₖ, Eₖ
    # return Hₖ_rows
end

𝜈(m) = ceil(m/2)

phases = range(0, π, length=50) # values of the adiabatic phase in (S32)
n_bands = 30
ee, EE = compute_bands_exact(;n_bands, phases, s, l, gₗ, Vₗ, λₗ, λₛ, ω)

fig1 = plot();
for i in 1:n_bands
    plot!(phases, EE[i, :], fillrange=EE[n_bands+i, :], fillalpha=0.1, label="band $i", legend=:outerright);
end
xlabel!(L"\varphi_t = \varphi_x"*", rad"); ylabel!("Floquet quasi-energy"*L"\varepsilon_{k,m}")

fig2 = plot();
for i in 1:2n_bands
    plot!(phases, ee[i, :], fillrange=ee[2n_bands+i, :], fillalpha=0.1, label="band $i", legend=:outerright);
end
xlabel!(L"\varphi_x"*", rad"); ylabel!("Energy "*L"\epsilon_{k,m}")

ee = compute_bands_exact(;n_bands=10, phases=[0], s, l, gₗ, Vₗ, λₗ, λₛ, ω)
scatter(zeros(length(bands)), bands)