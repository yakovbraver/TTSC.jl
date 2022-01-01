# The driving script used to analyse Hamiltonians (9) and (13) from arXiv:2012.02783

using Plots, LaTeXStrings
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

include("SpacetimeHamiltonian.jl")

function 𝐻₀(p, x, params)
    p^2 + params[1]*sin(x)^2
end

function 𝐻(p, x, params, t)
    p^2 + params[1]*sin(x)^2 + p * params[2] * params[3] * sin(params[3]*t)
end

function 𝑉(p::Real, x::Real)
    p
end

V₀ = 4320.0; ω = 240.0; λ = 0.01;
s = 3
params = [V₀, λ, ω]
H = SpacetimeHamiltonian(𝐻₀, 𝐻, (π/2, π), (π, 3π/2), (2.5, 3.5), (4.5, 5.5), params, s)

function plot_actions(H::SpacetimeHamiltonian)
    figs = [plot() for _ in 1:4];
    x = range(0, 2π, length=50);
    figs[1] = plot(x, H.𝑈, xlabel=L"x", ylabel=L"U(x)=V_0\sin^{2}(x)", legend=false);
    # title!(L"V_0 = %$(round(V₀, sigdigits=4))");
    I = Dierckx.get_knots(H.𝐸)
    figs[2] = plot(I, H.𝐸(I), xlabel=L"I", ylabel=L"E", legend=false);
    figs[3] = plot(I, H.𝐸′, xlabel=L"I", ylabel=L"dE/dI", xlims=(I[1], I[end]), legend=false);
    figs[4] = plot(I, H.𝐸″, xlabel=L"I", ylabel=L"d^2E/dI^2", xlims=(I[1], I[end]), legend=false);
    lay = @layout [a{0.5w} grid(3,1)]
    plot(figs..., layout=lay)
    # savefig("H0.pdf")
end

plot_actions(H)

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation

Iₛ, M, coeffs = compute_parameters(H, Function[𝑉])

function plot_isoenergies(; M, λ, ω, pₛ, Iₛ, s)
    ϑ = range(0, 2π, length=50)
    I = vcat(0:2:30, 30.5:0.5:42)
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = (I[i]-Iₛ)^2/2M - λ*ω*abs(pₛ)*cos(s*ϑ[t])
    end
    levels = vcat(range(minimum(E), -20, length=20), range(-19, maximum(E), length=10))
    contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (13)", color=:viridis; levels)
    hline!([Iₛ], label=L"I_s = %$(round(Iₛ, sigdigits=4))", c=:white) |> display
    title!(L"\lambda = %$(round(λ, sigdigits=2))")
end

plot_isoenergies(; pₛ=coeffs[1], M, λ, ω, Iₛ, s)
savefig("lambda_0.025/isoenergies.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

fig = plot();
for i in [2:2:34; Iₛ; 36:39]
    I, Θ = compute_IΘ(H, i)
    scatter!(mod2pi.(Θ.+pi/2), I, xlabel=L"\theta", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
    # scatter!(Θ, I, xlabel=L"\Theta=\theta-\omega t/s", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
end
display(fig)

fig = plot();
title!(L"\lambda = 0.01")
savefig("lambda0.01.pdf")