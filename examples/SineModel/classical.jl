# A driving script for analysing classical Hamiltonian (B1) from https://doi.org/10.1103/PhysRevB.106.144301 (https://arxiv.org/abs/2206.14804)
using TTSC.Classical
import TTSC.SineModel as sm
using Plots, LaTeXStrings, ProgressMeter

plotlyjs()
theme(:dark, size=(800, 600))

function 𝐻₀(p, x, params)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2
end

function 𝐻(p, x, params, t)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2 +
    params[4]*sin(2x)^2*cos(2params[6]*t) + 
    params[5]*cos(2x)^2*cos(params[6]*t + pi/2)
end

function 𝑄ₛ(p::Real, x::Real)
    sin(2x)^2
end

function 𝑄ₗ(p::Real, x::Real)
    cos(2x)^2
end

l = 1
gₗ = -7640
Vₗ = -2
λₛ = 100; λₗ = 40; ω = 410
s = 2
params = [gₗ, l, Vₗ, λₛ, λₗ, ω]
H = ClassicalHamiltonian(𝐻₀, 𝐻, params, s, min_pos=(1.5, 2), max_pos=(2, 2.5))

function plot_actions(H::ClassicalHamiltonian)
    figs = [plot() for _ in 1:4];
    x = range(0, 2π, length=200);
    I = Dierckx.get_knots(H.𝐸)
    figs[1] = plot!(x, H.𝑈, xlabel=L"x", c=1, title=L"V(x)=\lambda e^{-\sigma l/2}\cosh\sigma x", label=L"V(x)")
    figs[2] = plot(I, H.𝐸(I), xlabel=L"I", ylabel=L"E");
    figs[3] = plot(I, H.𝐸′, xlabel=L"I", ylabel=L"dE/dI", legend=false);
    figs[4] = plot(I, H.𝐸″, xlabel=L"I", ylabel=L"d^2E/dI^2", legend=false);
    lay = @layout [a{0.5w} grid(3,1)]
    plot(figs..., layout=lay)
end

plot_actions(H)

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [2s, s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
Aₗ = abs(coeffs[2]); χₗ = angle(coeffs[2])

# Below is the analysis of quantised classical Hamiltonian (B11)

########## Periodic case

φₜ = range(0, 2π, length=61)
n_cells = s
gₗ = -2λₛ*Aₛ
Vₗ = 2λₗ*Aₗ

h = sm.UnperturbedHamiltonian(n_cells; M, gₗ, Vₗ, φₓ=-φₜ/2, maxband=2, isperiodic=true)
sm.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(φₜ, r, label=false)
end
plot!(xlabel=L"\phi_t", ylabel="Energy")

# Wannier centres
sm.compute_wanniers!(h; targetband=1, mixsubbands=false)
fig = plot();
for (i, ϕ) in enumerate(φₜ)
    scatter!(h.w.pos[:, i], fill(ϕ, 2n_cells); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\phi_t", cbtitle="Energy")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₜ))
p = Progress(length(φₜ), 1)
@gif for (i, ϕ) in enumerate(φₜ)
    U = @. -λₛ*Aₛ*cos(4x) + λₗ*Aₗ*cos(2x - ϕ) + H.𝐸(Iₛ) - ω/s*Iₛ
    plot(x, U, label=false, ylims=(-5610, -5575))
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-5610, -5575))
    for j in 1:size(w, 2)
        plot!(x, 4abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false)
    end
    next!(p)
end

########## Non-periodic case

h = sm.UnperturbedHamiltonian(n_cells; M, gₗ, Vₗ, φₓ=-φₜ/2, maxband=2, isperiodic=false)
sm.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(φₜ, r, label=false)
end
plot!(xlabel=L"\phi_t", ylabel="Energy")

# Wannier centres
sm.compute_wanniers!(h; targetband=1, mixsubbands=false)

fig = plot();
for (i, ϕ) in enumerate(φₜ)
    scatter!(h.w.pos[:, i], fill(ϕ, size(h.w.pos, 1)); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi_t", cbtitle="Energy")

x = range(0, n_cells*π, length=50n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₜ))
p = Progress(length(φₜ), 1)
@gif for (i, ϕ) in enumerate(φₜ)
    U = @. -λₛ*Aₛ*cos(4x) + λₗ*Aₗ*cos(2x - ϕ) + H.𝐸(Iₛ) - ω/s*Iₛ
    plot(x, U, label=false, ylims=(-5610, -5575))
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-5610, -5575))
    for j in 1:size(w, 2)
        plot!(x, 4abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false)
    end
    next!(p)
end