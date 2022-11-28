using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

includet("DeltaModel.jl")
import .DeltaModel

### Unperturbed Hamiltonian

n_cells = 2
a = 6; λ = 500; U = 3
φₓ = range(0, 2π, length=61)
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)

# dispersion

function plot_dispersion(ε::AbstractVector; φ::Real, uh::DeltaModel.UnperturbedHamiltonian)
    cos_kL = [DeltaModel.cos_ka_tm(E; φ, uh) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(ka)", ylims=(-4, 4), ticks=:native, xlims=(0, ε[end]),
         title=L"U=%$U, a=%$a, \lambda=%$λ, \varphi=%$(round(φ, digits=3))", titlepos=:left, label=false)
    hline!([-1, 1], c=:white, label=false)
    vline!(((1:30) .* pi ./ (uh.a/3)).^2, c=2, label=L"(\frac{\pi n}{a/3})^2", lw=0.5) # analytical energy for a single well of width `a/3`
end

ε = range(U, 2000, length=2000)
plot_dispersion(ε; φ=φₓ[1], uh=h)
xlims!(340, 360)
savefig("dispersion.pdf")

import IntervalRootFinding as iroots
using IntervalArithmetic: (..)

f(E) = DeltaModel.cos_ka(E; φ=0, uh=h)
bounds = (50, 2000)
rts = iroots.roots(f, bounds[1]..bounds[2])
z = [rts[i].interval.lo for i in eachindex(rts)]
sort!(z)
scatter!(z, zeros(length(z)))

DeltaModel.diagonalise!(h, length(z), bounds)

### Floquet Hamiltonian

λₛ = 100; λₗ = 40; ω = 200
s = 2
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:space)

DeltaModel.diagonalise!(H)

skipbands = 2 # number of spatial band that have been skipped by the choice if `bounds` above
fig = plot();
for ik in axes(H.E, 2)
    for m in axes(H.E, 1)
        plot!(φₓ, H.E[m, ik, :] .- ω/s*skipbands)
    end
end
title!(L"N=%$n_cells, a=%$a, U=%$U, \lambda=%$(h.λ)")

# Ordered quasienergy spectrum
E_ordered = DeltaModel.order_floquet_levels(H)
fig = plot();
for ik in axes(E_ordered, 2)
    for m in axes(E_ordered, 1)
        plot!(φₓ, E_ordered[m, ik, :], label="band $(H.ν[m])", c=H.ν[m])
    end
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Maps of Floquet modes
x = range(0, n_cells*pi, length=50n_cells)
Ωt = range(0, 2π, length=40s)
iϕ = 1
whichstates = 1:4
u = DeltaModel.make_eigenfunctions(H, x, Ωt, [iϕ], whichstates) .|> abs2
figs = [plot() for _ in eachindex(whichstates)]
for (f, n) in enumerate(whichstates)
    figs[f] = heatmap(x, Ωt, u[:, :, n, iϕ]', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode $n")
end
plot(figs...)