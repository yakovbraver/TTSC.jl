using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

includet("DeltaModel.jl")
import .DeltaModel

### Unperturbed Hamiltonian

n_cells = 2
a = 4; λ = 10000; U = 1
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

ε = range(U, 5100, step=0.05)
plot_dispersion(ε; φ=φₓ[1], uh=h)
xlims!(340, 360)
savefig("dispersion.pdf")

import IntervalRootFinding as iroots
using IntervalArithmetic: (..)

f(E) = DeltaModel.cos_ka(E; φ=0, uh=h)
bounds = (45, 5100)
rts = iroots.roots(f, bounds[1]..bounds[2])
z = [rts[i].interval.lo for i in eachindex(rts)]
sort!(z)
scatter!(z, zeros(length(z)))

DeltaModel.diagonalise!(h, length(z), bounds)

# spectrum

fig = plot();
for m in axes(h.E, 2)
    for ik in axes(h.E, 1)
        band = ceil(Int, m/3)
        plot!(φₓ, h.E[ik, m, :], c=band, label="$band")
    end
end
title!(L"N=%$n_cells, a=%$a, U=%$U, \lambda=%$(h.λ)")

### Floquet Hamiltonian

# λₛ = 500; λₗ = 100; ω = 485
# 496
# 514
λₛ = 100; λₗ = 8; ω = 494
s = 2
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:both)

DeltaModel.diagonalise!(H)

skipbands = 2 # number of spatial band that have been skipped by the choice if `bounds` above
fig = plot();
for ik in axes(H.E, 2)
    for m in axes(H.E, 1)
        plot!(φₓ, H.E[m, ik, :] .- ω/s*skipbands, c=m)
    end
end
title!(L"N=%$n_cells, a=%$a, U=%$U, \lambda=%$(h.λ), \omega=%$ω")

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
n_x = 50
Ωt = range(0, 2π, length=40s)
iϕ = 1
whichsubbands = 1:3
x, u = DeltaModel.make_eigenfunctions(H, n_x, Ωt, [iϕ], whichsubbands)
u_real = abs2.(u) |> real
figs = [plot() for _ in 1:length(whichsubbands)*n_cells]
for (f, n) in enumerate(whichsubbands)
    for ik in 1:n_cells
        figs[n_cells*(f-1)+ik] = heatmap(x, Ωt, u_real[:, :, ik, n, iϕ]', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode n=$n, ik=$ik")
    end
end
plot(figs...)