using Plots, LaTeXStrings, ProgressBars
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

### Calculate isoenergies

include("effective.jl")

function plot_isoenergies(H::EffectiveHamiltonian)
    ϑ = range(0, 2π, length=50)
    I = range(H.I0-10, H.I0+10, length=50)
    E = energies_in_phasespace(H, ϑ, I)
    # calculate fixed energy limits so that the z-scale is the same regardless of the phases
    Elims = H.h0 - H.Ω*H.I0 .+ (-H.λL - H.λS, H.λL + H.λS + (I[end]-H.I0)^2/2H.M)
    fig = contourf(ϑ, I, E', xlabel=L"\theta"*", rad", ylabel=L"I", cbartitle="Energy", levels=range(Elims[1], Elims[2], length=20), color=:viridis)
    # fig = surface(ϑ, I, E', xlabel=L"\theta"*", rad", ylabel=L"I", z="Energy", color=:viridis)
    title!(L"M = %$(round(H.M, sigdigits=2)),"*
           L"\lambda_L = %$(round(H.λL, sigdigits=2)), \chi_L = %$(round(H.χL, sigdigits=2)),"*
           L"\lambda_S = %$(round(H.λS, sigdigits=2)), \chi_S = %$(round(H.χS, sigdigits=2))")
    fig
end

e = plot_isoenergies(EffectiveHamiltonian(h0=0.0, Ω=0.0, I0=1.0, M=2, λL=10, χL=pi, λS=5, χS=0.0, s=2))

ϕ = range(0, 2π, length=50)
H = EffectiveHamiltonian(h0=0.0, Ω=0.0, I0=0.0, M=0.5, λL=10, χL=0, λS=10, χS=0.0, s=2)
@gif for i in ProgressBar(ϕ)
    H.χL = i
    plot_isoenergies(H)
end

### Calculate bands

phases = range(0, 2π, length=40)
n_bands = 4
M = 0.04; λL = 50; λS = 50;
bands = get_bands(; n_bands, phases, M, λL, λS)

fig = plot()
for i in 1:n_bands
    plot!(fig, phases, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i")
end
xlabel!(L"\varphi_t"*", rad")
ylabel!("Energy")

### Calculate the spatial Hamiltonian as a function of the action, and compute the associated derivatives

import Polynomials

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

include("spatial.jl")

h0 = SpatialHamiltonian(l=2, g=4.0, V_L=1.0)
I, E = H_of_I(h0)

ws = 2 # window half-size for taking derivatives
E′, E″ = d_and_d²(I, E)

figs = [plot() for _ in 1:3];
figs[1] = plot(I, E, xlabel=L"I", ylabel=L"E", legend=false);
figs[2] = plot(I[ws+1:end-ws], E′[ws+1:end-ws], xlabel=L"I", ylabel=L"dE/dI", xlims=(I[1], I[end]), legend=false);
figs[3] = plot(I[ws+1:end-ws], E″[ws+1:end-ws], xlabel=L"I", ylabel=L"d^2E/dI^2", xlims=(I[1], I[end]), legend=false);
plot(figs..., layout=grid(3,1))