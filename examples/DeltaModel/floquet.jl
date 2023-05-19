# A driving script for analysing Floquet Hamiltonian (S1) from https://arxiv.org/abs/2305.07668
import TTSC.DeltaModel as dm
using Plots, LaTeXStrings, ProgressMeter

theme(:dark, size=(800, 500))

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

### Unperturbed Hamiltonian

n_cells = 2
a = 4; λ = 2000; U = 7
φₓ = range(0, 2π, length=31)
h = dm.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)

# dispersion

function plot_dispersion(ε::AbstractVector; φ::Real, uh::dm.UnperturbedHamiltonian)
    cos_kL = [dm.cos_ka_tm(E; φ, uh) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(ka)", ylims=(-4, 4), ticks=:native, xlims=(0, ε[end]),
         title=L"U=%$U, a=%$a, \lambda=%$λ, \varphi=%$(round(φ, digits=3))", titlepos=:left, label=false)
    hline!([-1, 1], c=:white, label=false)
    vline!(((1:30) .* pi ./ (uh.a/3)).^2, c=2, label=L"(\frac{\pi n}{a/3})^2", lw=0.5) # analytical energy for a single well of width `a/3`
end

ε = range(U, 6200, step=0.05)
plot_dispersion(ε; φ=φₓ[6], uh=h)

dm.diagonalise!(h; bounds=(300, 10000))

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

λₛ = 20; λₗ = 10; ω = 676.8
s = 2
pumptype = :spacetime
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)

# Quasienergy spectrum
skipbands = 7 # number of spatial bands that have been skipped by the choice if `bounds` above
fig = plot();
n_levels = size(h.E, 2)
for m in 1:n_levels
    for ik in 1:n_cells
        plot!(fig, φₓ, H.E[m, ik, :] .- ω/s*skipbands, label="sb $m (b $(H.ν[m]+skipbands)), k $ik", c=H.ν[m], xlabel=L"\varphi_x", ylabel="Quasienergy", ticks=:native)
    end
end
plot!(xlabel=L"\varphi_x", title=L"\omega=%$ω, \lambda_S=%$λₛ, \lambda_L=%$λₗ")

# Maps of Floquet modes
n_x = 50
Ωt = range(0, 2π, length=40s)
iφ = 1
whichsubbands = 11:12
x, u = dm.make_eigenfunctions(H, n_x, Ωt, [iφ], whichsubbands)
figs = [plot() for _ in 1:length(whichsubbands)*n_cells]
for n in eachindex(whichsubbands)
    for ik in 1:n_cells
        figs[n_cells*(n-1)+ik] = heatmap(x, Ωt, abs2.(u[:, :, ik, n, 1])', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode n=$n, ik=$ik")
    end
end
plot(figs...)

# Wannier centres
targetsubbands = 11:12
dm.compute_wanniers!(H; targetsubbands, slide_time=true)
fig = plot();
for (i, φ) in enumerate(φₓ)
    scatter!(H.w.pos[:, i], fill(φ, length(targetsubbands)*n_cells); label=false, markerstrokewidth=0, c=1)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x")

# Maps of Wannier functions
iφ = 1
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, [iφ])
figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
for f in eachindex(figs)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Wannier $f")
end
plot(figs..., plot_title="φ = $(round(φₓ[iφ], sigdigits=3))")

# Pumping animation

YELLOW = colorant"rgb(219, 173, 106)"
RED    = colorant"rgb(237, 71, 74)"
GREEN  = colorant"rgb(132, 221, 99)"
BLUE   = colorant"rgb(54, 201, 198)"
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)

function shadecells!(fig)
    fillalpha = 0.3; alpha = 0;
    top = [2, 1.5, 1, 0.5]
    bot = [1.5, 1, 0.5, 0]
    colours = [GREEN, RED, YELLOW, BLUE]
    for i in 0:4
        for j in 1:4
            plot!(fig, [1/6 + i*1/3, 1/6 + (i+1)*1/3], [top[j], top[j]]; fillrange=bot[j], fillalpha, alpha, c=colours[j])
        end
        colours .= circshift(colours, 2)
    end
    for j in 1:4
        plot!(fig, [0, 1/6], [top[j], top[j]]; fillrange=bot[j], fillalpha, alpha, c=colours[j])
        plot!(fig, [11/6, 2], [top[j], top[j]]; fillrange=bot[j], fillalpha, alpha, c=colours[j])
    end
    vline!([i/3 for i in 1:5], c=:black, lw=0.3)
end

"Return proper order of states depending on pumping type. (Hardcoded for `length(φₓ) == 31`.)"
function get_order(iφ::Integer; pumptype::Symbol)
    if pumptype == :time
        if iφ < 5
            order = [1, 2, 3, 4]
        elseif iφ < 20
            order = [2, 1, 4, 3]
        else
            order = [1, 2, 3, 4]
        end
    elseif pumptype == :space
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
    else
        if iφ < 4
            order = [1, 2, 3, 4]
        elseif iφ < 6
            order = [2, 1, 4, 3]
        elseif iφ < 16
            order = [3, 4, 1, 2]
        elseif iφ < 20
            order = [4, 3, 2, 1]
        elseif iφ < 26
            order = [3, 4, 1, 2]
        else
            order = [4, 3, 2, 1]
        end
    end
    return order
end

"Plot Wanniers at phase `iφ`. Pass values `iφ` > `length(φₓ)` to go beyound one period."
function plot_wanniers(w, iφ)
    figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
    phase = iφ > length(φₓ) ? iφ - length(φₓ) : iφ
    for (f, o) in enumerate(get_order(phase; pumptype))
        if iφ > length(φₓ)
            f = (f == 1 ? 4 : f == 4 ? 1 : f == 2 ? 3 : 2)
        end
        figs[f] = heatmap(x/a, Ωt/π, abs2.(w[:, :, o, phase]'), xlabel=L"x/a", ylabel=L"\Omega t/\pi", c=CMAP, title="Wannier $f", clims=(0, 3), xlims=(0, n_cells), ylims=(0, 2))
        shadecells!(figs[f])
    end
    plot(figs..., plot_title=L"\varphi_x=\varphi_t = %$(round(φₓ[phase]+2pi*(iφ > length(φₓ)), sigdigits=3))", plot_titlelocation=:left, legend=false)
end

gr()
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, eachindex(φₓ))

p = Progress(2length(φₓ), 1)
@gif for iφ in 1:2length(φₓ)
    plot_wanniers(w, iφ)
    next!(p)
end