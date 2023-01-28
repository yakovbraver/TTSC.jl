using Plots, LaTeXStrings, ProgressMeter
import IntervalRootFinding as iroots
using IntervalArithmetic: (..)
using LinearAlgebra: diagind

plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

includet("DeltaModel.jl")
import .DeltaModel

# Unperturbed Hamiltonian

n_cells = 1
a = 4; λ = 10000; U = 1
φₓ = range(0, 2π, length=31)
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)

f(E) = DeltaModel.cos_ka(E; φ=0, uh=h)
bounds = (45, 5100)
rts = iroots.roots(f, bounds[1]..bounds[2])
z = [rts[i].interval.lo for i in eachindex(rts)]

DeltaModel.diagonalise!(h, length(z), bounds)

# Floquet Hamiltonian

λₛ = 40; λₗ = 20; ω = 499.5
s = 2
pumptype = :time
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H)

skipbands = 2 # number of spatial bands that have been skipped by the choice if `bounds` above
fig = plot();
n_levels = size(h.E, 2)
for m in 1:n_levels
    for ik in 1:n_cells
        plot!(φₓ, H.E[m, ik, :] .- ω/s*skipbands, label="sb $m (b $(H.ν[m])), k $ik", c=H.ν[m], xlabel=L"\varphi_x", ylabel="Quasienergy", ticks=:native)
    end
end
plot!(xlabel=L"\varphi_x", title=L"\omega=%$ω, \lambda_S=%$λₛ, \lambda_L=%$λₗ")
savefig("ordered-spectrum.pdf")

### Floquet TB

# Construct Wanniers for the first temporal band
startsubband = 55
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)

DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)

# Maps of Wannier functions
iφ = 1
DeltaModel.optimise_wanniers!(H; whichstates=[(2, 3), (6, 7), (10, 11)], iφ)

CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true);
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [iφ])
figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
for f in eachindex(figs)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=CMAP, title="Wannier $f")
end
plot(figs..., layout=(length(figs)÷4, 4), plot_title=L"\lambda=%$λ, \omega=%$ω, %$iφ") |> display
savefig("$λ-$λₛ-$λₗ-wanniers.png")

DeltaModel.swap_wanniers!(H.w, 1, 2, eachindex(φₓ))
DeltaModel.swap_wanniers!(H.w, 6, 7, eachindex(φₓ))
DeltaModel.swap_wanniers!(H.w, 5, 6, eachindex(φₓ))
DeltaModel.swap_wanniers!(H.w, 9, 10, eachindex(φₓ))

# Construct the TB Floquet Hamiltonian

pumptype = :time
Htb = DeltaModel.TBFloquetHamiltonian(H; startsubband, pumptype)

# plot the Hamiltonian matrix
M = copy(Htb.H[:, :, 1])
M[diagind(M)] .= 0
M = abs.(M)
heatmap(1:size(M, 1), 1:size(M, 1), M; yaxis=:flip, c=CMAP, xticks=1:size(M, 1), yticks=1:size(M, 1), title=L"\lambda_S=%$λₛ, \lambda_L=%$λₗ, \omega=%$ω, \lambda=%$λ")
fontsize = 6
annotate!([(i, j, Plots.text(round(M[i, j], sigdigits=2), fontsize, RGB(0, 1, 0), :center)) for j in axes(M, 2) for i in axes(M, 1)])

plot(φₓ, real(Htb.H[1, 1, :]), label=L"H_{11}")
plot!(φₓ, real(Htb.H[2, 2, :]), label=L"H_{22}", title="λₛ = $λₛ, λₗ = $λₗ")
plot!(xlabel=L"\varphi_t", ylabel="(Quasi) Energy", title="Temporal on-site (quasi) energies, constant basis")

plot(φₓ, abs.(Htb.H[1, 7, :]), label=L"|H_{17}|")
plot!(φₓ, abs.(Htb.H[1, 11, :]), label=L"|H_{1,11}|")
plot!(xlabel=L"\varphi_x", ylabel="Coupling strength", title="Spatial couplings")

plot(φₓ, abs.(Htb.H[1, 2, :]), label=L"|H_{12}|")
plot!(φₓ, abs.(Htb.H[1, 4, :]), label=L"|H_{14}|")
plot!(xlabel=L"\varphi_t", ylabel="Coupling strength", title="Temporal couplings, constant basis")

# plot the path in parameter space (relevant only for the temporal pumping)
scatter(real(Htb.H[1, 1, :]) - real(Htb.H[2, 2, :]), abs.(Htb.H[1, 2, :]) - abs.(Htb.H[1, 4, :]), zcolor=φₓ, markerstrokewidth=0, markersize=6, c=:viridis,
        colorbar_title=L"\varphi_t", label="", xlabel=L"\Delta = H_{11} - H_{22}", ylabel=L"\delta = |H_{12}| - |H_{14}|",
        title=L"\lambda_S=%$λₛ, \lambda_L=%$λₗ, \omega=%$ω, \lambda=%$λ"*", constant basis", xlims=(-1, Inf))
savefig("$λₛ-$λₗ-bad-temporal-topology.pdf")

DeltaModel.diagonalise!(Htb)

skipbands = 2
fig = plot();
for r in eachrow(Htb.E)
    plot!(φₓ, r .- ω/s*skipbands)
end
plot!(xlabel=L"\varphi_t", ylabel="Quasienergy", title="TB", legend=false)

"Return proper order of states depending on pumping type."
function get_order_tb(iφ::Integer; pumptype::Symbol)
    if pumptype == :time
        if iφ < 17
            order = [1, 2]
        else
            order = [2, 1]
        end
    else
        if iφ < 6
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

function plot_tb_wanniers(w, iφ)
    figs = [plot() for _ in 1:length(targetlevels)]
    for (f, o) in enumerate(get_order_tb(iφ; pumptype))
        figs[f] = heatmap(1:3n_cells, 1:n_t, w[:, :, o, iφ]; clims, c=CMAP, xticks=1:3n_cells, yticks=1:n_t,
                          xlabel="spatial site", ylabel="temporal site", title="Wannier $f")
    end
    plot(figs..., plot_title=L"\varphi_x = %$(round(φₓ[iφ], sigdigits=3))")
end

targetlevels = [5, 6]
DeltaModel.compute_wanniers!(Htb; targetlevels)
whichphases = eachindex(φₓ)
w = abs2.(DeltaModel.make_wannierfunctions(Htb, whichphases))

clims = extrema(w)
n_t = 4 # number of temporal sites

iφ = 1
plot_tb_wanniers(w, iφ)

p = Progress(length(whichphases), 1)
@gif for iφ in whichphases
    plot_tb_wanniers(w, iφ)    
    next!(p)
end

# TB eigenstates
clims = extrema(abs2, Htb.c[:, :, iφ])
figs = [plot() for _ in 1:12]
for i in 1:12
    figs[i] = heatmap(1:3n_cells, 1:n_t, abs2.(reshape(Htb.c[:, i, iφ], (n_t, 3n_cells))); clims, c=CMAP, xticks=1:3n_cells, yticks=1:n_t,
                            xlabel="spatial site", ylabel="temporal site", title="Eigenstate $i")
end
plot(figs..., plot_title=L"\lambda_S=%$λₛ, \lambda_L=%$λₗ, \omega=%$ω, \lambda=%$λ")