using Plots, LaTeXStrings, ProgressBars
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

### Calculate the spatial Hamiltonian as a function of the action, and compute the associated derivatives

include("SpacetimeHamiltonian.jl")

function 𝐻₀(x, p, params)
    p[1]^2 + params[1]*sin(x[1])^2
end

function d𝑥╱d𝑡!(dx, x, p, params, t)
    dx[1] = 2p[1] + params[2] * params[3] * sin(params[3]*t)
end

function d𝑝╱d𝑡!(dp, x, p, params, t)
    dp[1] = -params[1] * sin(2x[1])
end

function 𝑉(x::Real, p::Real)
    p
end

V₀ = 4320.0; ω = 240.0; λ = 0.01;
s = 3 # freely chosen parameters
params = [V₀, λ, ω]
H = SpacetimeHamiltonian(𝐻₀, (π/2, π), (π, 3π/2), (2.5, 3.5), (4.5, 5.5), d𝑥╱d𝑡!, d𝑝╱d𝑡!, params, s)

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

### Set main parameters
Iₛ, M, coeffs = compute_pₛ(H, Function[𝑉])
@code_warntype compute_pₛ(H, Function[𝑉])
### Calculate isoenergies

function plot_isoenergies(ϑ::AbstractVector, I::AbstractVector; M, λ, ω, pₛ, Iₛ, s)
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = (I[i]-Iₛ)^2/2M - λ*ω*abs(pₛ)*cos(s*ϑ[t])
    end
    levels = vcat(range(minimum(E), -20, length=20), range(-19, maximum(E), length=10))
    contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (13)", color=:viridis; levels)
    hline!([Iₛ], label=L"I_s = %$(round(Iₛ, sigdigits=4))", c=:white) |> display
    title!(L"\lambda = %$(round(λ, sigdigits=2))")
end

ϑ = range(0, 2π, length=50)
I = vcat(0:2:30, 30.5:0.5:42)

plot_isoenergies(ϑ, I; M, λ, ω, pₛ, Iₛ, s)
savefig("lambda_0.025/isoenergies.pdf")

### Calculate evolutions of Hamiltonian (11)

fig = plot();
for I in [2:2:34; Iₛ; 36:39]
    I, Θ = compute_IΘ(H, I)
    scatter!(Θ, I, xlabel=L"\theta", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
    # scatter!(Θ, I, xlabel=L"\Theta=\theta-\omega t/s", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
end
display(fig)

fig = plot();
title!(L"\lambda = 0.01")
savefig("lambda0.01.pdf")