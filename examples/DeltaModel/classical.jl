# A driving script for the classical analysis of Hamiltonian (2) from https://arxiv.org/abs/2305.07668
using TTSC.Classical
using Plots, LaTeXStrings, ProgressMeter

plotlyjs()
theme(:dark, size=(800, 600))

function 𝐻₀(p, x, params)
    σ, l, λ = params
    return p^2 + λ * exp(-σ*l/2) * cosh(σ*x)
end

function 𝐻(p, x, params, t)
    σ, l, λ, λₛ, λₗ, ω = params
    return p^2 + λ * exp(-σ*l/2) * cosh(σ*x) + 
           λₛ * 𝑄ₛ(p, x) * cos(2ω*t) + 
           λₗ * 𝑄ₗ(p, x) * cos(ω*t)
end

function 𝑄ₛ(p::Real, x::Real)
    cos(12π*x/a)
end

function 𝑄ₗ(p::Real, x::Real)
    cos(6π*x/a)
end

a = 4.0
l = a/3
σ = 100 / l
λ = 1000
λₛ = 20; λₗ = 10; ω = 676.8
s = 2
params = [σ, l, λ, λₛ, λₗ, ω]

H = ClassicalHamiltonian(𝐻₀, 𝐻, params, s, left_tp=(-1.06l/2, 0.0), right_tp=(0.0, 1.06l/2))

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [2s, s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
Aₗ = abs(coeffs[2]); χₗ = angle(coeffs[2])

import Dierckx
function plot_actions(H::ClassicalHamiltonian)
    figs = [plot() for _ in 1:4];
    x = range(-1.2l/2, 1.2l/2, length=200);
    I = Dierckx.get_knots(H.𝐸)
    figs[1] = vline([-l/2, l/2], c=:white, label=L"x = \pm l", legendposition=(0.5, 0.8))
    plot!(x, H.𝑈, xlabel=L"x", c=1, title=L"V(x)=\lambda e^{-\sigma l/2}\cosh\sigma x", label=L"V(x)", ylims=(-100, 5000))
    figs[2] = plot(I, H.𝐸(I), xlabel=L"I", ylabel=L"E", label="numerical", legendposition=:topleft, ylims=(-10, 6000));
    plot!(I, (π*I/l).^2, label="exact");
    figs[3] = plot(I, H.𝐸′, xlabel=L"I", ylabel=L"dE/dI", legend=false);
    plot!(I, 2(π/l)^2*I, label="exact");
    figs[4] = plot(I, H.𝐸″, xlabel=L"I", ylabel=L"d^2E/dI^2", legend=false);
    hline!([2(π/l)^2], label="exact");
    lay = @layout [a{0.5w} grid(3,1)]
    plot(figs..., layout=lay)
end

plot_actions(H)
savefig("h_0-parameters.pdf")

### Check accuracy and stability

# sample initial conditions
x₀ = -1.04l/2
p₀ = 0.0

# calculate `n_T` periods of unperturbed motion to check accuracy
import OrdinaryDiffEq as DiffEq
using DiffEqPhysics: HamiltonianProblem

T =  l/sqrt(𝐻₀(p₀, x₀, params)) # analytical period 
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

I_min = 30; I_max = 40
ϑ = range(0, 2π, length=100)
I = range(I_min, I_max, length=50)
E = Matrix{Float64}(undef, length(ϑ), length(I))
h₀ = H.𝐸(Iₛ) - ω/s*Iₛ
for i in eachindex(I), t in eachindex(ϑ)
    E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] - χₛ) + λₗ*Aₗ*cos(s*ϑ[t] - χₗ)
end
figa = contour(ϑ ./ π, I, E', xlabel=L"\Theta/\pi", ylabel=L"I", minorgrid=true, c=[:white], colorbar=false)
savefig(figa, "secular.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

function point_to_angle(p, x, E, T)
    if x > 0 && p > 0
        t = (x + l/2) / (2√E) - T / 4
    elseif p < 0
        t = (x - l/2) / (-2√E) + T / 4
    else
        t = (x + l/2) / (2√E) + 3T / 4
    end
    return 2π * t / T
end

figb = plot();
for (I_min, χ₀) in zip([30, 30, 30], [0, 1, -1])
    for i in I_min:0.2:I_max
        display(i)
        I, Θ = compute_IΘ(H, i; n_T=100, χ₀, point_to_angle)
        scatter!(Θ ./ π, I, xlabel=L"\Theta/\pi", ylabel=L"I", markerstrokewidth=0, label=false, minorgrid=true, markersize=2)
    end    
end
xlims!(0, 2)
savefig(figb, "exact.pdf")

# Quantisation

import TTSC.SineModel as sm

φₜ = range(0, 2π, length=61)
n_cells = s
gₗ = 2λₛ*Aₛ
Vₗ = -2λₗ*Aₗ
M = l^2 / 2π^2 # analytical result

h = sm.UnperturbedHamiltonian(n_cells; M, gₗ=gₗ, Vₗ=Vₗ, φₓ=-φₜ/2, maxband=4, isperiodic=true)
sm.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(φₜ, r, label=false)
end
plot!(xlabel=L"\varphi_t", ylabel="Energy")
savefig("qc-spectrum.pdf")

# Wannier centres

sm.compute_wanniers!(h; targetband=1, mixsubbands=false)
fig = plot();
for (i, φ) in enumerate(φₜ)
    scatter!(h.w.pos[:, i], fill(φ, 2n_cells); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"\Theta", ylabel=L"\varphi_t", cbtitle="Energy")
savefig("qc-centres.pdf")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₜ))
p = Progress(length(φₜ), 1)
@gif for (i, φ) in enumerate(φₜ)
    U = @. gₗ*cos(4x) + Vₗ*cos(2x - φ) + H.𝐸(Iₛ) - ω/s*Iₛ
    plot(x, U, label=false, ylims=(h.w.E[1, 1]-10, h.w.E[3, 1]+10))
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0, colorbar=false)
    for j in axes(w, 2)
        plot!(x, 4abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false, xlabel=L"\Theta")
    end
    next!(p)
end