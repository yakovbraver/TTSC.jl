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
φₓ = range(0, 2π, length=31)
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
        plot!(φₓ, h.E[ik, m, :], c=band, label="$band", xlabel=L"\varphi_x", ylabel=L"\varepsilon", legend=false)
    end
end
display(fig)

### Floquet Hamiltonian

λₛ = 40; λₗ = 20; ω = 499.5
s = 2
pumptype = :space
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H)

# Quasienergy spectrum
skipbands = 2 # number of spatial bands that have been skipped by the choice if `bounds` above
fig = plot();
n_levels = size(h.E, 2)
for m in 1:n_levels
    for ik in 1:n_cells
        plot!(fig, φₓ, H.E[m, ik, :] .- ω/s*skipbands, label="sb $m (b $(H.ν[m])), k $ik", c=H.ν[m], xlabel=L"\varphi_x", ylabel="Quasienergy", ticks=:native)
    end
end
plot!(xlabel=L"\varphi_x", title=L"\omega=%$ω, \lambda_S=%$λₛ, \lambda_L=%$λₗ")

# Maps of Floquet modes
n_x = 50
Ωt = range(0, 2π, length=40s)
iφ = 1
whichsubbands = range(61, length=3)
x, u = DeltaModel.make_eigenfunctions(H, n_x, Ωt, [iφ], whichsubbands)
figs = [plot() for _ in 1:length(whichsubbands)*n_cells]
for n in eachindex(whichsubbands)
    for ik in 1:n_cells
        figs[n_cells*(n-1)+ik] = heatmap(x, Ωt, abs2.(u[:, :, ik, n, 1])', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode n=$n, ik=$ik")
    end
end
plot(figs...)

# Wannier centres
targetsubbands = [60, 63]
DeltaModel.compute_wanniers!(H; targetsubbands, slide_time=true)
fig = plot();
for (i, φ) in enumerate(φₓ)
    scatter!(H.w.pos[:, i], fill(φ, length(targetsubbands)*n_cells); label=false, markerstrokewidth=0, c=1)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x")

# Maps of Wannier functions
iφ = 1
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [iφ])
figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
for f in eachindex(figs)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Wannier $f")
end
plot(figs..., plot_title="φ = $(round(φₓ[iφ], sigdigits=3))")

# Pumping animation

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
        if iφ < 6
            order = [1, 2, 3, 4]
        elseif iφ < 8
            order = [4, 3, 2, 1]
        elseif iφ < 17
            order = [3, 4, 1, 2]
        elseif iφ < 26
            order = [4, 3, 2, 1]
        else
            order = [3, 4, 1, 2]
        end
    end
    return order
end
pyplot()
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, eachindex(φₓ))
p = Progress(length(φₓ), 1)
@gif for iφ in eachindex(φₓ)
    figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
    for (f, o) in enumerate(get_order(iφ; pumptype))
        figs[f] = heatmap(x, Ωt/π, abs2.(w[:, :, o, iφ]'), xlabel=L"x", ylabel=L"\Omega t/\pi", c=CMAP, title="Wannier $f", clims=(0, 3), xlims=(0, a*n_cells), ylims=(0, 2))
        shadecells!(figs[f])
    end
    plot(figs..., plot_title=L"\varphi_x=\varphi_t = %$(round(φₓ[iφ], sigdigits=3))")
    next!(p)
end