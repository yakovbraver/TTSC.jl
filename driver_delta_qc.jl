using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
theme(:dark, size=(800, 600))

include("SpacetimeHamiltonian.jl")

function 𝐻₀(p, x, params)
    σ, l, λ = params
    return p^2 + λ*exp(-σ*l/2)cosh(σ*x)
end

function 𝐻(p, x, params, t)
    σ, l, λ, λₛ, λₗ, ω = params
    return p^2 + λ*exp(-σ*l/2)cosh(σ*x) + 
           λₛ * sin(2x)^2 * cos(2ω*t) + 
           λₗ * cos(2x)^2 * cos(ω*t)
end

function 𝑄ₛ(p::Real, x::Real)
    cos(2x)^2
end

function 𝑄ₗ(p::Real, x::Real)
    cos(2x)^2
end

function f0!(dp, p, x, params, t)
    σ, l, λ = params
    dp[1] = - λ*σ*exp(-σ*l/2)sinh(σ*x[1])
end

function f1!(dp, p, x, params, t)
    σ, l, λ, λₛ, λₗ, ω = params
    dp[1] = - λ*σ*exp(-σ*l/2)sinh(σ*x[1]) + 
              2λₛ * sin(4x[1]) * cos(2ω*t) + 
              2λₗ * sin(4x[1]) * cos(ω*t + pi/2)
end

function f2!(dx, p, x, params, t)
    dx[1] = 2p[1]
end

σ = 150.0
a = 2.0
params = [σ, a]
x₀ = -1.054
p₀ = 0.0
tspan = (0, 2)
H_problem = DynamicalODEProblem(f0!, f2!, [p₀], [x₀], tspan, params)
sol = DiffEq.solve(H_problem, DiffEq.KahanLi8(), dt=2e-4)
sol = DiffEq.solve(H_problem, DiffEq.McAte5(), dt=2e-4)

plot(sol)
plot(sol, vars=2)
plot(sol.t, sol[1, :])

a = 4.0
l = a/3
σ = 100 / l
λ = 500
λₛ = 100; λₗ = 40; ω = 300
s = 2
params = [σ, l, λ, λₛ, λₗ, ω]

H_problem = DynamicalODEProblem(f1!, f2!, [p₀], [x₀], tspan, params)
sol = DiffEq.solve(H_problem, DiffEq.KahanLi8(), dt=2e-4)
plot(sol, vars=2)

H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, left_tp=(-1.05l/2, 0.0), right_tp=(0.0, 1.05l/2))

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [2s, s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
Aₗ = abs(coeffs[2]); χₗ = angle(coeffs[2])

function plot_actions(H::SpacetimeHamiltonian)
    figs = [plot() for _ in 1:4];
    x = range(-1.2l/2, 1.2l/2, length=200);
    I = Dierckx.get_knots(H.𝐸)
    figs[1] = vline([-l/2, l/2], c=:white, label=L"x = \pm l", legendposition=(0.5, 0.8))
    plot!(x, H.𝑈, xlabel=L"x", c=1, title=L"V(x)=\lambda e^{-\sigma l/2}\cosh\sigma x", label=L"V(x)", ylims=(-100, 5000))
    figs[2] = plot(I, H.𝐸(I), xlabel=L"I", ylabel=L"E", label="numerical", legendposition=:topleft, ylims=(-10, 5000));
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

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation
theme(:default, size=(800, 600))

I_min = 16; I_max = 28
ϑ = range(0, 2π, length=100)
I = range(I_min, I_max, length=50)
E = Matrix{Float64}(undef, length(ϑ), length(I))
h₀ = H.𝐸(Iₛ) - ω/s*Iₛ
for i in eachindex(I), t in eachindex(ϑ)
    E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] - χₛ) + λₗ*Aₗ*cos(s*ϑ[t] - χₗ + pi/2)
end
figa = contour(ϑ ./ π, I, E', xlabel=L"\Theta/\pi", ylabel=L"I", minorgrid=true, c=[:white], colorbar=false)
ylims!(figa, (17.5, 30))
savefig(figa, "secular.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space
plotlyjs()
pyplot()
figb = plot();
for (I_min, χ₀) in zip([16, 22.5, 23.5], [0, 0.5, -0.5])
    for i in I_min:0.5:I_max
        display(i)
        I, Θ = compute_IΘ(H, i; n_T=100, χ₀)
        scatter!(Θ ./ π, I, xlabel=L"\Theta/\pi", markerstrokewidth=0, label=false, minorgrid=true, markersize=2)
    end    
end
xlims!(0, 2)
savefig(figb, "exact.pdf")

plot(figa, figb)

# Quantisation

include("bandsolvers.jl")

import .Bandsolvers

φₜ = range(0, 2π, length=61)
n_cells = s
χₛ
χₗ
gₗ = 2λₛ*Aₛ
Vₗ = 2λₗ*Aₗ

h = Bandsolvers.UnperturbedHamiltonian(n_cells; M=M/9, gₗ=gₗ, Vₗ=Vₗ, φₓ=-φₜ/2, maxband=2, isperiodic=true)
Bandsolvers.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(φₜ, r, label=false)
end
plot!(xlabel=L"\phi_t", ylabel="Energy")
savefig("qc-spectrum.pdf")

# Wannier centres
plotlyjs()
Bandsolvers.compute_wanniers!(h; targetband=1)
fig = plot();
for (i, ϕ) in enumerate(φₜ)
    scatter!(h.w.pos[:, i], fill(ϕ, 2n_cells); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\phi_t", cbtitle="Energy")
savefig("qc-centres.pdf")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
_, w = Bandsolvers.make_wannierfunctions(h, x, 1:length(φₜ))
p = Progress(length(φₜ), 1)
@gif for (i, ϕ) in enumerate(φₜ)
    U = @. λₛ*Aₛ*cos(4x) + λₗ*Aₗ*cos(2x - ϕ) + H.𝐸(Iₛ) - ω/s*Iₛ
    plot(x, U, label=false, ylims=(h.w.E[1, 1]-10, h.w.E[3, 1]+10))
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:size(w, 2)
        plot!(x, 4abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false)
    end
    next!(p)
end