using Plots, LaTeXStrings, ProgressBars
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

### Calculate the spatial Hamiltonian as a function of the action, and compute the associated derivatives

import Polynomials

include("spatial.jl")

"""
Return the first and second derivatives of `y(x)` as a tuple `(y′, y″)`. Use window half-size `ws`.
The values in `y′` and `y″` will correspond to the vector `x`; first and last `ws` values will contain zeros.
"""
function d_and_d²(x, y; ws=2)
    N = length(y)
    y′ = zeros(N)
	y″ = zeros(N)
	for i in 1+ws:N-ws
		f = Polynomials.fit(x[i-ws:i+ws], y[i-ws:i+ws], 2)
		y′[i] = Polynomials.derivative(f, 1)(x[i])
		y″[i] = Polynomials.derivative(f, 2)(x[i])
	end
	y′, y″
end

l = 2; g = 40.0; V_L = 0.0
h0 = SpatialHamiltonian(; l, g, V_L)
I, E = H_of_I(h0)

ws = 2 # window half-size for taking derivatives
E′, E″ = d_and_d²(I, E)

figs = [plot() for _ in 1:4];
x = range(0, π, length=30)
figs[1] = plot(x, x -> 𝑈(h0, x), xlabel=L"x", ylabel=L"U(x)=\tilde{g}_\ell\cos^{2\ell}(2x)+V_L\cos^{2}(x)", legend=false);
title!(L"\ell = %$l, g = %$(round(g, sigdigits=2)), V_L = %$(round(V_L, sigdigits=2))");
figs[2] = plot(I, E, xlabel=L"I", ylabel=L"E", legend=false);
figs[3] = plot(I[ws+1:end-ws], E′[ws+1:end-ws], xlabel=L"I", ylabel=L"dE/dI", xlims=(I[1], I[end]), legend=false);
figs[4] = plot(I[ws+1:end-ws], E″[ws+1:end-ws], xlabel=L"I", ylabel=L"d^2E/dI^2", xlims=(I[1], I[end]), legend=false);
lay = @layout [a{0.5w} grid(3,1)]
plot(figs..., layout=lay)
savefig("h0-(S2).pdf")

### Calculate isoenergies

using Interpolations: LinearInterpolation
include("effective.jl")

function plot_isoenergies(H::EffectiveHamiltonian, ϑ::AbstractVector, I::AbstractVector)
    E = energies_in_phasespace(H, ϑ, I)
    # Calculate fixed energy limits so that the z-scale is the same regardless of the phases χL and χS.
    # For the potential λL⋅cos(𝑠ϑ + χL) + λS⋅cos(2𝑠ϑ + χS),
    # the lowest possible value is -λL-λS (attainable at some ϑ when |χL - χS| = π/2 + π𝑛),
    # and the highest is λL+λS (attainable at some ϑ when |χL - χS| = 2π𝑛)
    if H.M < 0
        Elims = H.h0 - H.Ω*H.I0 .+ (-H.λL - H.λS + (I[end]-H.I0)^2/2H.M,    H.λL + H.λS)
    else
        Elims = H.h0 - H.Ω*H.I0 .+ (-H.λL - H.λS,    H.λL + H.λS + (I[end]-H.I0)^2/2H.M)
    end
    contourf(ϑ, I, E', xlabel=L"\theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (S17)", levels=range(Elims[1], Elims[2], length=20), color=:viridis)
    hline!([H.I0], label=L"I_0 = %$(round(H.I0, sigdigits=2))", c=:black)
    # surface(ϑ, I, E', xlabel=L"\theta"*", rad", ylabel=L"I", zlabel="Energy", zlims=(Elims[1], Elims[2]), color=:viridis)
    title!(L"M = %$(round(H.M, sigdigits=2)),"*
           L"\lambda_L = %$(round(H.λL, sigdigits=2)), \chi_L = %$(round(H.χL, sigdigits=2)),"*
           L"\lambda_S = %$(round(H.λS, sigdigits=2)), \chi_S = %$(round(H.χS, sigdigits=2))")
end

E_of_I  = LinearInterpolation(I, E)
E′_of_I = LinearInterpolation(I, E′)
E″_of_I = LinearInterpolation(I, E″)

I0 = 2.53 # the working point of choice
E0 = E_of_I(I0) # energy at the working point
Ω = E′_of_I(I0) # frequency at the working point
M = 1 / E″_of_I(I0) # "mass" at the working point
λL = 8; χL = 0; λS = 5; χS = 0.0; s = 2 # other freely chosen parameters for (S17)
Heff = EffectiveHamiltonian(h0=E0; Ω, I0, M, λL, χL, λS, χS, s)

ϑ = range(0, 2π, length=50)
I = range(0, 2I0, length=50)
plot_isoenergies(Heff, ϑ, I)

ϕ = range(0, 2π, length=30)
anim = @animate for i in ProgressBar(ϕ)
    Heff.χL = i
    plot_isoenergies(Heff, ϑ, I)
end
gif(anim, "isoenergies.gif", fps = 10)

### Calculate bands

ϕₜ = range(0, 2π, length=40) # values of the adiabatic phase in (S32)
n_bands = 4
bands = get_bands(; n_bands, phases=ϕₜ, M=abs(M), λL, λS)

plot();
for i in 1:n_bands
    plot!(ϕₜ, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i");
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Energy")
savefig("bands.pdf")