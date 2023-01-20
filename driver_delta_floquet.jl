using Plots, LaTeXStrings, ProgressMeter
import IntervalRootFinding as iroots
using IntervalArithmetic: (..)

plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

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
        plot!(φₓ, h.E[ik, m, :], c=band, #=label="$band",=# xlabel=L"\varphi_x", ylabel=L"\varepsilon", legend=false)
    end
end
title!(L"N=%$n_cells, a=%$a, U=%$U, \lambda=%$(h.λ)")
ylims!(2934, 2938)
savefig("spatial-spectrum-zoom.pdf")

### Floquet Hamiltonian

λₛ = 10; λₗ = 5; ω = 494
s = 2
pumptype = :time
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H)

skipbands = 2 # number of spatial bands that have been skipped by the choice if `bounds` above
fig = plot();
for ik in axes(H.E, 2)
    for m in axes(H.E, 1)
        plot!(φₓ, H.E[m, ik, :] .- ω/s*skipbands, c=m)
    end
end
plot!(title=L"\omega=%$ω, \lambda_S=%$λₛ, \lambda_L=%$λₗ", ticks=:native)

# Ordered quasienergy spectrum
E_ordered = DeltaModel.order_floquet_levels(H)
fig = plot();
for ik in axes(E_ordered, 2)
    for m in axes(E_ordered, 1)
        plot!(φₓ, E_ordered[m, ik, :], label="band $(H.ν[m])", c=H.ν[m], xlabel=L"\varphi_x", ylabel="Quasienergy", legend=false, ticks=:native)
    end
end
plot!(xlabel=L"\varphi_t", title=L"\omega=%$ω, \lambda_S=%$λₛ, \lambda_L=%$λₗ")
savefig("floquet-space.pdf")

# Maps of Floquet modes
n_x = 50
Ωt = range(0, 2π, length=40s)
iφ = 1
whichsubbands = 4:6
x, u = DeltaModel.make_eigenfunctions(H, n_x, Ωt, [iφ], whichsubbands)
u_real = abs2.(u)
figs = [plot() for _ in 1:length(whichsubbands)*n_cells]
for n in eachindex(whichsubbands)
    for ik in 1:n_cells
        figs[n_cells*(n-1)+ik] = heatmap(x, Ωt, u_real[:, :, ik, n, 1]', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode n=$n, ik=$ik")
    end
end
plot(figs...)

# Wannier centres
targetsubbands = [3, 6]
DeltaModel.compute_wanniers!(H; targetsubbands)
fig = plot();
for (i, φ) in enumerate(φₓ)
    scatter!(H.w.pos[:, i], fill(φ, length(targetsubbands)*n_cells); label=false, markerstrokewidth=0, c=1)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x")
savefig("centres-space.pdf")

# Maps of Wannier functions
iφ = 12
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [iφ])
figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
for f in eachindex(figs)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Wannier $f")
end
plot(figs..., plot_title="φ = $(round(φₓ[iφ], sigdigits=3))")

# Pumping animation

theme(:default, size=(800, 500))
GREEN = colorant"rgb(132, 221, 99)"
RED   = colorant"rgb(237, 71, 74)"
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)

function shadecells!(fig)
    fillalpha = 0.3; alpha = 0;
    x = a/6
    for i in 0:2:6n_cells-1
        plot!(fig, [i*x, (i+1)x], [2, 2]; fillrange=1.7, fillalpha, alpha, c=GREEN)
        plot!(fig, [i*x, (i+1)x], [0.7, 0.7]; fillrange=0, fillalpha, alpha, c=GREEN)
        plot!(fig, [(i+1)x, (i+2)x], [1.7, 1.7]; fillrange=0.7, fillalpha, alpha, c=GREEN, legend=false)
        vline!([i*x], c=RED)
    end
end

"Return proper order of states depending on pumping type."
function get_order(iφ::Integer; pumptype::Symbol)
    if pumptype == :time
        order = [1, 2, 3, 4]
    else
        if iφ < 12
            order = [1, 2, 3, 4]
        elseif iφ < 31
            order = [4, 3, 2, 1]
        elseif iφ < 51
            order = [3, 4, 1, 2]
        else
            order = [4, 3, 2, 1]
        end
    end
    order
end

n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, 1:length(φₓ))
p = Progress(length(φₓ), 1)
@gif for iφ in eachindex(φₓ)
    # figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
    figs = [plot() for _ in 1:4]
    for (f, o) in enumerate(get_order(iφ; pumptype))
        figs[f] = heatmap(x, Ωt/π, abs2.(w[:, :, o, iφ]'), xlabel=L"x", ylabel=L"\Omega t/\pi", c=CMAP, title="Wannier $f", clims=(0, 3), xlims=(0, a*n_cells), ylims=(0, 2))
        shadecells!(figs[f])
    end
    plot(figs..., plot_title=L"\varphi_x=\varphi_t = %$(round(φₓ[iφ], sigdigits=3))")
    next!(p)
end

### Floquet TB

## Construct Wanniers for the first temporal band

λₛ = 10; λₗ = 5; ω = 494
s = 2
pumptype = :time

n_cells = 1
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)
f(E) = DeltaModel.cos_ka(E; φ=0, uh=h)
bounds = (45, 5100)
rts = iroots.roots(f, bounds[1]..bounds[2])
DeltaModel.diagonalise!(h, length(rts), bounds)
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:time)
DeltaModel.diagonalise!(H)

targetsubbands = 1:12
n_subbands = length(targetsubbands)

DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.3π)
DeltaModel.optimise_wanniers!(H, [(2, 3), (6, 7), (10, 11)], iφ=1)

# Maps of Wannier functions
iφ = 1
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [1])
figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
# figs = [plot() for _ in 1:4]
for f in eachindex(figs)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=CMAP, title="Wannier $f")
end
plot(figs...)
savefig("all-wanniers-12.png")

DeltaModel.swap_wanniers!(H.w, 3, 1)
DeltaModel.swap_wanniers!(H.w, 2, 3)
DeltaModel.swap_wanniers!(H.w, 5, 6)
DeltaModel.swap_wanniers!(H.w, 11, 9)
DeltaModel.swap_wanniers!(H.w, 10, 11)

# Construct the TB Floquet Hamiltonian

Htb = DeltaModel.TBFloquetHamiltonian(H; N=n_cells, isperiodic=true, targetband=1, pumptype=:space)

# plot the Hamiltonian matrix
using LinearAlgebra: diagind
M = copy(Htb.H[:, :, 1])
M[diagind(M)] .= 0
heatmap(1:size(M, 1), 1:size(M, 1), abs.(M); yaxis=:flip, c=CMAP, xticks=1:size(M, 1), yticks=1:size(M, 1))

# Htb = DeltaModel.TBFloquetHamiltonian(H, isperiodic=true)

DeltaModel.diagonalise!(Htb)

skipbands = 2
fig = plot();
for r in eachrow(Htb.E)
    plot!(φₓ, r .- ω/s*skipbands)
end
plot!(xlabel=L"\varphi_t", ylabel="Quasienergy", title="TB", legend=false, ylims=(-2752, -2728))

skipbands = 2 # number of spatial bands that have been skipped by the choice if `bounds` above
fig = plot();
for ik in axes(H.E, 2)
    for m in axes(H.E, 1)
        plot!(φₓ, H.E[m, ik, :] .- ω/s*skipbands, c=1)
    end
end
plot!(xlabel=L"\varphi_t", ylabel="Quasienergy", legend=false, ylims=(-2752, -2728))

targetsubbands = [5, 6]
DeltaModel.compute_wanniers!(Htb; targetsubbands)
whichphases = 1:length(φₓ)
w = abs2.(DeltaModel.make_wannierfunctions(Htb, whichphases))

"Return proper order of states depending on pumping type."
function get_order_tb(iφ::Integer; pumptype::Symbol)
    if pumptype == :time
        if iφ < 17
            order = [1, 2]
        elseif iφ < 20
            order = [2, 1]
        elseif iφ < 28
            order = [1, 2]
        else
            order = [2, 1]
        end
    else
        if iφ < 5
            order = [1, 2]
        elseif iφ < 16
            order = [2, 1]
        elseif iφ < 27
            order = [1, 2]
        else
            order = [2, 1]
        end
    end
    order
end

clims = extrema(w)
n_t = size(Htb.H, 1) ÷ 3n_cells # number of temporal sites
p = Progress(length(whichphases), 1)
@gif for iφ in whichphases
    figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
    for (f, o) in enumerate(get_order_tb(iφ; pumptype=:space))
        figs[f] = heatmap(1:3n_cells, 1:n_t, w[:, :, o, iφ]; clims, c=CMAP, xticks=1:3n_cells, yticks=1:n_t,
                          xlabel="spatial site", ylabel="temporal site", title="Wannier $f")
    end
    # plot(figs..., plot_title=L"\varphi_x=\varphi_t = %$(round(φₓ[i], sigdigits=3))")
    plot(figs..., plot_title=L"\varphi_t = %$(round(φₓ[iφ], sigdigits=3))")
    next!(p)
end