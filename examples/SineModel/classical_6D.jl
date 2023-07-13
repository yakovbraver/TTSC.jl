# A driving script for analysing classical Hamiltonian (S2) from https://doi.org/10.1103/PhysRevB.103.L100301 (https://arxiv.org/abs/2012.02783)
using TTSC.Classical
using Plots, LaTeXStrings

plotlyjs()
theme(:dark, size=(800, 500))

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
# plot(range(0, 2π, length=200), x -> 𝐻₀(0, x, params))
H = ClassicalHamiltonian(𝐻₀, 𝐻, params, s, min_pos=(3.0, 3.2), max_pos=(1.4, 1.6))

import Dierckx
function plot_actions(H::ClassicalHamiltonian)
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
end

plot_actions(H)
savefig("h_0-parameters.pdf")

### Check accuracy and stability

# sample initial conditions
p₀ = 0.0; x₀ = 2.0;

# Calculate analytical period (ZP-2.2.3.18).
# In our case the equation of motion is
#     𝑞″ + 2𝑉₀sin(2𝑞) = 0,
# and we choose the initial conditions 𝑞(0) = 𝑥₀, 𝑞′(0) = 0.
# We let 2𝑞 = 𝑦 to get the standard form
#     𝑦″ + 4𝑉₀sin(𝑦) = 0
# with the initial condition 𝑦(0) = 2𝑥₀. The period is then
#     𝑇 = 4 / √(4𝑉₀) 𝐾(𝑚²), where 𝑚 = sin(𝑦(0)/2) = sin(𝑥₀)
import OrdinaryDiffEq as DiffEq
using DiffEqPhysics: HamiltonianProblem
import Elliptic
m = sin(x₀)
T = 2 / √V₀ * Elliptic.K(m^2)

# calculate `n_T` periods of unperturbed motion to check accuracy
n_T = 100
tspan = (0, n_T*T)
H_problem = HamiltonianProblem(𝐻₀, p₀, x₀, tspan, params)
sol = DiffEq.solve(H_problem, DiffEq.McAte3(); dt=2e-4)
plot(sol)
vline!([i*T for i = 1:n_T], c=:white)

# calculate `n_T` periods of perturbed motion to check stability
H_problem = HamiltonianProblem(𝐻, p₀, x₀, tspan, params)
sol = DiffEq.solve(H_problem, DiffEq.McAte3(); dt=2e-4)
plot(sol)

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation

Iₛ, M, coeffs = compute_parameters(H, Function[𝑉], [s])

function plot_isoenergies(; M, λ, ω, pₛ, Iₛ, s)
    ϑ = range(0, 2π, length=50)
    I_max = last(Dierckx.get_knots(H.𝐸))
    I = [0:2:30; range(30.5, I_max, length=20)]
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = (I[i]-Iₛ)^2/2M - λ*ω*abs(pₛ)*cos(s*ϑ[t])
    end
    levels = [range(minimum(E), -20, length=20); range(-19, maximum(E), length=10)]
    contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (13)", color=:viridis; levels)
    hline!([Iₛ], label=L"I_s = %$(round(Iₛ, sigdigits=4))", c=:white)
    title!(L"\lambda = %$(round(λ, sigdigits=2))")
end

pₛ = abs(coeffs[1])
plot_isoenergies(; pₛ, M, λ, ω, Iₛ, s)
savefig("secular-isoenergies.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

fig = plot();
for i in 25:0.5:40
    display(i)
    I, Θ = compute_IΘ(H, i, χ₀=1, n_T=200)
    scatter!(mod2pi.(Θ.+pi/2), I, xlabel=L"\theta", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false, markersize=2)
end
display(fig)
savefig("exact-isoenergies.pdf")