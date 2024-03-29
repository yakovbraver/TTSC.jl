# A driving script for analysing Floquet Hamiltonian (A1) from https://doi.org/10.1103/PhysRevB.106.144301 (https://arxiv.org/abs/2206.14804)
import TTSC.SineModel as sm
using Plots, LaTeXStrings, ProgressMeter

plotlyjs()
theme(:dark, size=(800, 600))

l = 1
gₗ = -7640
Vₗ = -2
λₛ = 100; λₗ = 40; ω = 410
s = 2

########## Periodic case

φₓ = [range(0, pi/4-0.1, length=10); range(pi/4-0.01, pi/4+0.01, length=10);
      range(pi/4+0.1, 3pi/4-0.1, length=20); range(3pi/4-0.01, 3pi/4+0.01, length=10);
      range(3pi/4+0.1, pi, length=10)]
n_cells = 2

h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=30, isperiodic=true)
sm.diagonalise!(h)

# unperturbed Hamiltonian spectrum
fig = plot();
plot!(range(0, π, length=200), x -> gₗ*cos(2x)^2 + Vₗ*cos(x)^2, lw=2, c=:white, label=false) # spatial potential
for r in eachrow(h.E)
    plot!(φₓ, r, label=false)
end
plot!(xlabel=L"\phi_x", ylabel="Energy")

H = sm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:both, minband=1)
sm.diagonalise!(H)

# Floquet quasienergy spectrum
fig = plot();
for r in eachrow(H.E)
    plot!(φₓ, r, label=false)
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Ordered quasienergy spectrum
E_ordered = sm.order_floquet_levels(H)
fig = plot();
for (i, r) in enumerate(eachrow(E_ordered))
    m = i + H.minlevel - 1
    plot!(φₓ, r, label="band $(H.ν[m]), level $m", c=H.ν[m])
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Maps of Floquet modes
x = range(0, n_cells*pi, length=50n_cells)
Ωt = range(0, 2π, length=40s)
iϕ = 15
whichstates = 1:4
u = sm.make_eigenfunctions(H, x, Ωt, [iϕ], whichstates) .|> abs2
figs = [plot() for _ in eachindex(whichstates)]
for (f, n) in enumerate(whichstates)
    figs[f] = heatmap(x, Ωt, u[:, :, n, 1]', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode $n")
end
plot(figs...)

# Wannier centres
targetlevels = [1, 2, 5, 6]
sm.compute_wanniers!(H; targetlevels)
fig = plot();
for (i, ϕ) in enumerate(φₓ)
    scatter!(H.uh.w.pos[:, i], fill(ϕ, length(targetlevels)); label=false, markerstrokewidth=0, c=1)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\phi_x")

# Maps of Wannier functions
_, w = sm.make_wannierfunctions(H, x, Ωt, [iϕ])
figs = [plot() for _ in eachindex(targetlevels)]
for f in eachindex(targetlevels)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Wannier $f")
end
plot(figs...)

### TB Floquet Hamiltonian

# Construct the Wanniers
targetlevels = 1:4*2n_cells
sm.compute_wanniers!(H; targetlevels)

# Plot the Wanniers
x = range(0, n_cells*pi, length=50n_cells)
Ωt = range(0, 2π, length=40s)
iϕ = 1
_, w = sm.make_wannierfunctions(H, x, Ωt, [iϕ])
figs = [plot() for _ in eachindex(targetlevels)]
for f in eachindex(targetlevels)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Wannier $f")
end
plot(figs...)

# Construct and diagonalise TB Hamiltonian
Htb = sm.TBFloquetHamiltonian(H; targetband=1, pumptype=:space)
sm.diagonalise!(Htb)

# Quasienergy spectrum
fig = plot();
for r in eachrow(Htb.E)
    plot!(φₓ, r, label=false)
end
plot!(xlabel=L"\varphi_x", ylabel="Quasienergy", title="TB")

########## Non-periodic case

h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=30, isperiodic=false)
sm.diagonalise!(h)

# unperturbed Hamiltonian spectrum
fig = plot();
plot!(range(0, π, length=200), x -> gₗ*cos(2x)^2 + Vₗ*cos(x)^2, lw=2, c=:white, label=false) # spatial potential
for r in eachrow(h.E)
    plot!(φₓ, r, label=false)
end
plot!(xlabel=L"\phi_x", ylabel="Energy")

H = sm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:space, minband=1)
sm.diagonalise!(H)

# Floquet quasienergy spectrum
fig = plot();
for r in eachrow(H.E)
    plot!(φₓ, r, label=false)
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Ordered quasienergy spectrum
E_ordered = sm.order_floquet_levels(H)
fig = plot();
for (i, r) in enumerate(eachrow(E_ordered))
    m = i + H.minlevel - 1
    plot!(φₓ, r, label="band $(H.ν[m]), level $m", c=H.ν[m])
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Maps of Floquet modes
x = range(0, n_cells*pi, length=50n_cells)
Ωt = range(0, 2π, length=40s)
iϕ = 1
whichstates = 1:3
u = sm.make_eigenfunctions(H, x, Ωt, [iϕ], whichstates) .|> abs2
figs = [plot() for _ in eachindex(whichstates)]
for (f, n) in enumerate(whichstates)
    figs[f] = heatmap(x, Ωt, u[:, :, n, 1]', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode $n")
end
plot(figs...)