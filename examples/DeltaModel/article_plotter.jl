# A script for producing all the figures for https://doi.org/10.1103/physrevb.108.l020303 (https://arxiv.org/abs/2305.07668)
import TTSC.DeltaModel as dm
using Plots, Measures, LaTeXStrings
using LinearAlgebra

pythonplot()

"Set plotting defaults and initialise the canvas size with the given `width` and `height` (in cm)"
function set_defaults(;width, height)
    cm2px = 39
    theme(:default, size=(width*cm2px, height*cm2px))
    default(legendfontsize=8, tickfontsize=8, titlefontsize=9, labelfontsize=8, colorbar_titlefontsize=8, plot_titlefontsize=9, xwiden=false, frame=:box, annotationfontsize=10)
end

YELLOW = colorant"rgb(219, 173, 106)"
RED    = colorant"rgb(237, 71, 74)"
GREEN  = colorant"rgb(132, 221, 99)"
BLUE   = colorant"rgb(54, 201, 198)"
BLUE2  = colorant"rgb(56, 111, 164)"
BLACK  = colorant"rgb(13, 19, 33)"
ORANGE = colorant"rgb(240, 161, 37)"
GREEN2 = colorant"rgb(69, 169, 58)"
BLUE3  = colorant"rgb(82, 152, 197)"
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

##########
########## FIG S2
##########

set_defaults(width=2*8.6, height=7.5)

### (a) Floquet spectrum

n_cells = 2
a = 4; λ = 2000; U = 7
n_φₓ = 401 # use many phases to resolve the anticrossings in the Floquet spectrum
φₓ = range(0, 2pi, n_φₓ)
h = dm.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)
dm.diagonalise!(h; bounds=(300, 10000))
skipbands = 7 # number of spatial bands that have been skipped by the choice of `bounds` above
λₛ = 20; λₗ = 10; ω = 676.8
s = 2
pumptype = :spacetime
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=true)

ind = [44:91; 178:201]
dm.swap_eigenstates!(H, 69, 74, [1], [ind; n_φₓ+1 .- ind])
ind = [21:24; 111:132; 144:157]
dm.swap_eigenstates!(H, 72, 74, [1], [ind; n_φₓ+1 .- ind])
ind = [44:91; 178:201]
dm.swap_eigenstates!(H, 71, 73, [1], [ind; n_φₓ+1 .- ind])
ind = 183:201
dm.swap_eigenstates!(H, 64, 71, [1], [ind; n_φₓ+1 .- ind])
ind = [13:42; 48:87; 91:123; 146:174; 182:201]
dm.swap_eigenstates!(H, 64, 73, [1], [ind; n_φₓ+1 .- ind])
ind = [28:108; 161:201]
dm.swap_eigenstates!(H, 65, 72, [2], [ind; n_φₓ+1 .- ind])
ind = [19:117; 152:201]
dm.swap_eigenstates!(H, 68, 73, [2], [ind; n_φₓ+1 .- ind])

m1 = 64
figS2a = plot();
for m in m1:m1+11
    for ik in 1:n_cells
        c = (m == 66 || m == 75 ? ORANGE : BLACK)
        # label = "$m ($(H.ν[m])), $ik"
        plot!(φₓ ./ π, H.E[m, ik, :] .+ 8ω .- ω/s*skipbands; c, ylabel="Quasienergy", lw=0.3)
    end
end
plot!(xlabel=L"\varphi_x=\varphi_t\ (\pi"*" rad)", legend=false)
annotate!(figS2a, 0.5, 240, (L"\nu_2=-1", RED))
annotate!(figS2a, 0.5, 270, (L"\nu_2=+1", RED))

### (b)-(c) Wannier maps

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

function four_wanniers(w, φ_str, order)
    figs = [plot() for _ in 1:4]
    figs[1] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, order[1], 1])', xformatter=_->"", ylabel=L"\Omega t/\pi", title=L"|w_{1}(x,t)|^2", c=CMAP, cbar=false)
    figs[2] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, order[2], 1])', xformatter=_->"", yformatter=_->"", title=L"|w_{2}(x,t)|^2", c=CMAP, cbar=false)
    figs[3] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, order[3], 1])', xlabel=L"x/a", ylabel=L"\Omega t/\pi", title=L"|w_{3}(x,t)|^2", c=CMAP, cbar=false)
    figs[4] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, order[4], 1])', xlabel=L"x/a", yformatter=_->"", title=L"|w_{4}(x,t)|^2", c=CMAP, cbar=false)
    foreach(shadecells!, figs)
    return plot(figs..., layout=(2, 2), link=:both, plot_title=L"\varphi_x=\varphi_t="*φ_str, xaxis=((0, 2), 0:1:2), yaxis=((0, 2), 0:1:2), legend=false)
end

targetsubbands = [m1+2, m1+11]
dm.compute_wanniers!(H; targetsubbands)
iφ = 1
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, [iφ])
figb = four_wanniers(w, "0", 1:4)

iφ = 68
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, [iφ])
figc = four_wanniers(w, L"\pi/3", [3, 4, 2, 1])

plot(figS2a, figb, figc, layout=grid(1, 3, widths=[0.285, 0.35, 0.35]), plot_title="", dpi=600)
savefig("figS2.png")

##########
########## FIG S3 (Wanniers)
##########

set_defaults(width=2*8.6, height=7.5)

function wanniers12(w)
    figs = [plot() for _ in 1:12]
    figs[1] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, 1, 1])', xformatter=_->"", ylabel=L"\Omega t/\pi", title=L"|w_{1}(x,t)|^2", c=CMAP, cbar=false)
    figs[5] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, 5, 1])', xformatter=_->"", ylabel=L"\Omega t/\pi", title=L"|w_{5}(x,t)|^2", c=CMAP, cbar=false)
    figs[9] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, 21, 1])', xlabel=L"x/a", ylabel=L"\Omega t/\pi", title=L"|w_{21}(x,t)|^2", c=CMAP, cbar=false)
    for i in [2:4; 6:8]
        figs[i] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, i, 1])', xformatter=_->"", yformatter=_->"", title=L"|w_{%$i}(x,t)|^2", c=CMAP, cbar=false)
    end
    for i in 22:24
        figs[i-12] = heatmap(x ./ a, Ωt ./ π, abs2.(w[:, :, i, 1])', xlabel=L"x/a", yformatter=_->"", title=L"|w_{%$i}(x,t)|^2", c=CMAP, cbar=false)
    end
    foreach(shadecells!, figs)
    return plot(figs..., layout=(3, 4), link=:both, xaxis=((0, 2), 0:1:2), yaxis=((0, 2), 0:1:2), legend=false, dpi=600)
end

n_φₓ = 2 # we are only interested in 𝜑ₓ = 0
φₓ = range(0, 2pi, n_φₓ)
h = dm.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)
dm.diagonalise!(h; bounds=(300, 10000))
pumptype = :spacetime
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)

startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)

dm.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
dm.order_wanniers!(H; optimise=true)

iφ = 1
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, [iφ])
fig = wanniers12(w)
savefig("figS3.png")

##########
########## FIG S4 (TB Hamiltonian matrix)
##########

set_defaults(width=8.6, height=7.5)

# using `H` from Fig S3
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)

M = copy(Htb.H[:, :, 1])
M[diagind(M)] .= 0
M = abs.(M)
heatmap(1:size(M, 1), 1:size(M, 1), M; yaxis=:flip, c=CMAP, xticks=1:4:24, yticks=1:4:24,
        xlabel=L"\ell", ylabel=L"\ell'", cbar_title=L"|J_{\ell'\ell}|", rightmargin=-5mm)
savefig("figS4.pdf")

##########
########## FIG 1
##########

set_defaults(width=8.6, height=5)
# set_defaults(width=5, height=4) # for poster

function plot_potential(;U::Real, lift::Real=0, iφ::Integer)
    V = [U * cos(φₓ[iφ] + 2π*(site%3)/3) + lift for site in 0:3n_cells-1]
    x = range(0, n_cells, 3n_cells+1)
    plot(x, V, lw=2, c=BLACK, label=false, xlabel=L"x/a", ylabel="Quasienergy", lt=:steppost) # `lt=:steppost` isn't supported on `plotlyjs`
    vline!(x, c=BLACK, label=false, lw=2) # barriers
end

n_φₓ = 61
φₓ = range(0, 2pi, n_φₓ)
h = dm.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)
dm.diagonalise!(h; bounds=(300, 10000))

pumptype = :space
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)
H.E .+= 8ω - ω/s*skipbands
startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)
dm.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
dm.order_wanniers!(H; optimise=true)

iφ = 7
n_x = 500
Ωt = [7π/4]
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, [iφ])

colours = [GREEN2, RED, ORANGE, BLUE3]
fig = plot_potential(;U, lift=253, iφ)

for site in 0:3n_cells-1
    xmask = site*n_x+1:(site+1)*n_x+1 # select only x's from the curernt site so as not to draw the zero-valued tails
    for f in 1:2s
        i = 2s*site+f # wannier number (1 to 2s*3N)
        plot!(x[xmask] ./ a, abs2.(w[xmask, 1, i, 1]) .+ H.w.E[i, iφ], c=colours[f], legend=false, lw=0.3)
    end
end
fig
savefig("fig1.pdf")

# 2D Wanniers (does not appear in the paper)
plotlyjs()
set_defaults(width=15, height=10)
xmask = 1:n_x+1 # select only x's from the curernt site so as not to draw the zero-valued tails
figs = [plot() for _ in 1:2s*2s]
for i in 1:2s
    wi = abs2.(w[xmask, 1, i, 1])
    for j in 1:2s
        wj = abs2.(w[xmask, 1, j, 1])
        figs[(i-1)*2s + j] = heatmap(x[xmask], x[xmask], wi * wj', title="$i, $j", c=:viridis, cbar=false)
    end
end
plot(figs...)

##########
########## FIG 2
##########

set_defaults(width=2*8.6, height=10)

### (a) Separated spatial spectrum

# using `H` from Fig 1
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)

# spatial part
H_S = dm.separate_space(Htb; N=n_cells, periodic=true)
S_E = Matrix{Float64}(undef, size(H_S, 1), size(H_S, 3))
for iφ in axes(H_S, 3)
    S_E[:, iφ] = eigvals(Hermitian(H_S[:, :, iφ])) #.+ 8ω .- ω/s*skipbands
end
figa = plot();
for r in eachrow(S_E)
    plot!(φₓ ./ π, r, legend=false, xlabel=L"\varphi_x/\pi", ylabel="Energy", c=BLACK, title="Spectrum of "*L"\hat{H}_x", lw=0.3)
end
plot!(ylims=(245, 270))
annotate!(0.7, 250, (L"\nu_1=-1", RED))
annotate!(0.7, 257, (L"\nu_1=+2", RED))
annotate!(0.7, 267, (L"\nu_1=-1", RED))

### (b) Separated temporal spectrum

pumptype = :time
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)
dm.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
dm.order_wanniers!(H; optimise=false)
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)

T_E = Matrix{Float64}(undef, 4, size(Htb.H, 3))
for iφ in eachindex(φₓ)
    T_E[:, iφ] = eigvals(Htb.H[1:4, 1:4, iφ])
end
T_E .+= -T_E[1, 1] + (minimum(H.E[:, :, 1]) + 8ω - ω/s*skipbands) - S_E[1, 1]
figb = plot();
for r in eachrow(T_E)
    plot!(φₓ ./ π, r, legend=false, xlabel=L"\varphi_t/\pi", ylabel="Energy", c=BLACK, title="Spectrum of "*L"\hat{H}_t", lw=0.3)
end
plot!(ylims=(-12.5, 12.5))
annotate!(0.7, -8, (L"\nu_1=+1", RED));
annotate!(0.7,  4, (L"\nu_1=-1", RED))

### (c) H'_TB spectrum

sz = size(S_E, 1)*size(T_E, 1)
HTB = Matrix{Float64}(undef, sz, length(φₓ))
for i = axes(S_E, 1)
    for j = axes(T_E, 1)
        HTB[size(T_E, 1)*(i-1)+j, :] .= S_E[i, :] .+ T_E[j, :]
    end
end
figc = plot();
for r in eachrow(HTB)
    plot!(φₓ ./ π, r, legend=false, xlabel=L"\varphi_x=\varphi_t\ (\pi"*" rad)", ylabel="Energy", c=BLACK, title="Spectrum of "*L"\hat{{\cal H}}_{{\rm TB}}^{\prime}", lw=0.3)
end
annotate!(0.7, 240, (L"\nu_2=-1", RED))
annotate!(0.7, 270, (L"\nu_2=+1", RED))

## gaps in 8D sectrum if two copies of H'_TB are combined

# higher band
hi = maximum(S_E[6, :]) + maximum(T_E[4, :])
mid = minimum(S_E[5, :]) + minimum(T_E[3, :])
lo = maximum(S_E[4, :]) + maximum(T_E[4, :])
width = 2hi - 2mid
gap = 2mid - (lo + hi)
gap/width

# lower band
hi = minimum(S_E[3, :]) + minimum(T_E[2, :])
mid = maximum(S_E[2, :]) + maximum(T_E[2, :])
lo = minimum(S_E[1, :]) + minimum(T_E[1, :])
width = 2mid - 2lo
gap = (lo + hi) - 2mid
gap/width

# sum of two approximate Hamiltonians H'_TB (does not appear in the paper)
# sort the levels because `HTB` was construced by summing two spectra, and the levels are not sorted
for c in eachcol(HTB)
    sort!(c)
end
E4D = Matrix{Float64}(undef, 6, length(φₓ))
E4D[1, :] = 2HTB[1, :]
E4D[2, :] = 2HTB[4, :]
E4D[3, :] = HTB[1, :] .+ HTB[5, :]
E4D[4, :] = HTB[20, :] .+ HTB[24, :]
E4D[5, :] = 2HTB[21, :]
E4D[6, :] = 2HTB[24, :]

figd = plot();
plot!(φₓ ./ π, E4D[1, :], fillrange=E4D[2, :], fillalpha=0.3, alpha=0, c=BLACK);
plot!(φₓ ./ π, E4D[3, :], fillrange=E4D[4, :], fillalpha=0.3, alpha=0, c=BLACK);
plot!(φₓ ./ π, E4D[5, :], fillrange=E4D[6, :], xlabel=L"\varphi_x=\varphi_y=\varphi_{t_x}=\varphi_{t_y}\ (\pi"*" rad)",
      ylabel="Energy", c=BLACK, title="Spectrum of "*L"\hat{{\cal H}}_{8{\rm D}}", fillalpha=0.3, alpha=0, legend=false)
hspan!([maximum(E4D[2, :]), minimum(E4D[3, :])], c=RED, alpha=0.5);
hspan!([maximum(E4D[4, :]), minimum(E4D[5, :])], c=RED, alpha=0.5)

### (f) 4D/8D spectrum

pumptype = :spacetime
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)
dm.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
dm.order_wanniers!(H; optimise=false)
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)
dm.diagonalise!(Htb)
Htb.E .+= 8ω - ω/s*skipbands

E8D = Matrix{Float64}(undef, 6, length(φₓ))
E8D[1, :] = 2Htb.E[1, :]
E8D[2, :] = 2Htb.E[4, :]
E8D[3, :] = Htb.E[1, :] .+ Htb.E[5, :]
E8D[4, :] = Htb.E[20, :] .+ Htb.E[24, :]
E8D[5, :] = 2Htb.E[21, :]
E8D[6, :] = 2Htb.E[24, :]

figf = plot();
plot!(φₓ ./ π, E8D[1, :], fillrange=E8D[2, :], fillalpha=0.3, alpha=0, c=BLACK);
plot!(φₓ ./ π, E8D[3, :], fillrange=E8D[4, :], fillalpha=0.3, alpha=0, c=BLACK);
plot!(φₓ ./ π, E8D[5, :], fillrange=E8D[6, :], xlabel=L"\varphi_x=\varphi_y=\varphi_{t_x}=\varphi_{t_y}\ (\pi"*" rad)",
      ylabel="Energy", c=BLACK, title="Spectrum of "*L"\hat{{\cal H}}_{8{\rm D}}", fillalpha=0.3, alpha=0, legend=false)
hspan!([maximum(E8D[2, :]), minimum(E8D[3, :])], c=RED, alpha=0.5);
hspan!([maximum(E8D[4, :]), minimum(E8D[5, :])], c=RED, alpha=0.5)
annotate!(0.7, 482, (L"\nu_4=+1", RED))
annotate!(0.7, 542, (L"\nu_4=+1", RED))
plot!(figS2a, title="Spectrum of "*L"\hat{{\cal H}}_{{\rm TB}}", ylabel="Energy")
plot(figa, figb, figc, figS2a, figf, layout=grid(2, 4), plot_title="")
savefig("fig2.pdf")

### (e) 3D energy surfaces of H_TB

set_defaults(width=8.6, height=7.5)

pumptype = :spacetime
M = 2 * 3n_cells*s
S = Array{Float64, 3}(undef, M, length(φₓ), length(φₓ)) # Spectrum of H_TB in the format S[n, φₓ, φₜ]
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
using ProgressMeter
@showprogress for iφ in eachindex(φₓ)
    dm.diagonalise!(H; reorder=false, φₜ=fill(φₓ[iφ], length(φₓ)))
    S[1:end÷2, :, iφ] = H.E[1:M÷2, 1, :] # save all levels at ik = 1
    S[end÷2+1:end, :, iφ] = H.E[1:M÷2, 2, :] # save all levels at ik = 2
end
S.+= 8ω - ω/s*skipbands
fig2e = plot();
for n in axes(S, 1)
    surface!(φₓ ./π, φₓ ./π, S[n, :, :], xlabel=L"\varphi_t/\pi", ylabel=L"\varphi_x/\pi", cbar=false,
             camera=(130, 3), topmargin=0mm, bottommargin=0mm, leftmargin=0mm, dpi=600, c=:Paired_8)
end
fig2e
savefig("fig2e_new.png")

# higher gap of 8D system
maximum(S[[11:12; 23:24], :, :])
hi = maximum(S[[11:12; 23:24], :, :])
mid = minimum(S[[11:12; 23:24], :, :])
lo = maximum(S[[9:10; 21:22], :, :])
width = 2hi - 2mid
gap = 2mid - (lo + hi)
gap/width # 5%

# lower gap of 8D system
minimum(S[1:2, :, :])
minimum(S[[3:4; 15:16], :, :])
lo = minimum(S[[1:2; 13:14], :, :])
mid = maximum(S[[1:2; 13:14], :, :])
hi = minimum(S[[3:4; 15:16], :, :])
width = 2mid - 2lo
gap = (lo + hi) - 2mid
gap/width # 2%

##########
########## FIG S5
##########

set_defaults(width=8.6, height=5)

### (a) Spatial pumping

pumptype = :space
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)
startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)
dm.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
dm.order_wanniers!(H; optimise=true)
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)

H_S = dm.separate_space(Htb; N=n_cells, periodic=true)
figa = plot();
for (targetlevels, c) in zip([1:2, 3:4, 5:6], [BLACK, RED, BLUE2])
    pos = dm.compute_wanniers(H_S; targetlevels, isperiodic=true, n_s=3)
    for r in eachrow(pos)
        scatter!(r, φₓ ./ π; c, legend=false, markerstrokewidth=0, markersize=2)
    end
end
plot!(figa, xlabel=L"j", ylabel=L"\varphi_x/\pi", widen=false, xlims=(0, 6))

### (b) Temporal pumping

pumptype = :time
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)
dm.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
dm.order_wanniers!(H; optimise=true)
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)

figb = plot();
for (targetlevels, c) in zip([1:2, 3:4], [BLACK, BLUE2])
    pos = dm.compute_wanniers(Htb.H[1:4, 1:4, :]; targetlevels, isperiodic=true, n_s=2)
    for r in eachrow(pos)
        scatter!(r, φₓ ./ π; c, legend=false, markerstrokewidth=0, markersize=2)
    end
end
plot!(xlabel=L"\alpha", ylabel=L"\varphi_t/\pi", widen=false, xlims=(0, 4))

plot(figa, figb, layout=(1, 2))
savefig("figS5.pdf")