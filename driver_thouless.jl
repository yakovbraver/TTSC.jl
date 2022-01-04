using Plots, LaTeXStrings
using SpecialFunctions: gamma
using Combinatorics: factorial
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
# plot(range(0, 2π, length=200), x -> 𝐻₀(0, x, params))
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

using BandedMatrices: BandedMatrix, band
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
    H = BandedMatrix{ComplexF64}(undef, (2n_j + 1, 2n_j + 1), (2, 2))
    H[band(-2)] .= λₛAₛ
    H[band(2)]  .= λₛAₛ
    
    bands = Matrix{Float64}(undef, 2n_bands, length(phases))
    for k in [0, s÷2] # iterate over the centre of BZ and then the boundary
        H[band(0)] .= [(2j + k)^2 / M for j = -n_j:n_j]
        # `a` and `b` control where to place the eigenvalues depedning on `k`; see description of `bands`
        a = k*n_bands + 1 
        b = a+n_bands - 1
        for (i, ϕ) in enumerate(phases)
            H[band(-1)] .= λₗAₗ*cis(-ϕ)
            H[band(1)]  .= λₗAₗ*cis(ϕ)
            vals, _, _ = eigsolve(H, n_bands, :SR)
            bands[a:b, i] .= vals[1:n_bands]
        end
    end
    return bands / 4π # we do not include the factor 4π in the diagonalisation problem and restore it here
end

phases = range(0, 2π, length=40) # values of the adiabatic phase in (S32)
n_bands = 4
bands = compute_bands(; n_bands, phases, s, M=abs(M), λₗAₗ=λₗ*Aₗ, λₛAₛ=λₛ*Aₛ)

plot();
for i in 1:n_bands
    plot!(phases, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i");
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Energy")
savefig("bands.pdf")