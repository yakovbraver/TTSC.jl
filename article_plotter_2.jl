using Plots, Measures, LaTeXStrings
using LinearAlgebra

pythonplot()
plotlyjs()

"Set plotting defaults and initialise the canvas size with the given `width` and `height` (in cm)"
function set_defaults(;width, height)
    cm2px = 39
    theme(:default, size=(width*cm2px, height*cm2px))
    default(legendfontsize=8, tickfontsize=8, titlefontsize=9, labelfontsize=8, colorbar_titlefontsize=8, plot_titlefontsize=9, xwiden=false, frame=:box, annotationfontsize=10)
end

YELLOW = colorant"rgb(219, 173, 106)"
GREY   = colorant"rgb(105, 116, 124)"
RED    = colorant"rgb(237, 71, 74)"
GREEN  = colorant"rgb(132, 221, 99)"
BROWN  = colorant"rgb(78, 1, 16)"   # 224, 186, 215
BLUE   = colorant"rgb(54, 201, 198)"
BLUE2   = colorant"rgb(56, 111, 164)"
GREY2  = colorant"rgb(180, 180, 180)"
BLACK  = colorant"rgb(13, 19, 33)"
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

include("DeltaModel.jl")
import .DeltaModel

##########
########## FIG 1
##########

set_defaults(width=2*8.6, height=7.5)

### (a) Floquet spectrum

n_cells = 2
a = 4; λ = 2000; U = 7
φₓ = range(0, 2pi, length=31)
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)
DeltaModel.diagonalise!(h; bounds=(300, 10000))
skipbands = 7 # number of spatial bands that have been skipped by the choice if `bounds` above
λₛ = 20; λₗ = 10; ω = 676.8
s = 2
pumptype = :spacetime
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=true)

ind = [7; 15:17; 27]
DeltaModel.swap_eigenstates!(H, 69, 74, [1], [ind; 62 .- ind])
ind = 18:24
DeltaModel.swap_eigenstates!(H, 72, 74, [1], [ind; 62 .- ind])
ind = [18]
DeltaModel.swap_eigenstates!(H, 65, 72, [2], [ind; 62 .- ind])
ind = [1:2; 9:14; 20:22]
DeltaModel.swap_eigenstates!(H, 64, 73, [1], [ind; 62 .- ind])
ind = [8]
DeltaModel.swap_eigenstates!(H, 71, 73, [1], [ind; 62 .- ind])
ind = 9:14
DeltaModel.swap_eigenstates!(H, 64, 71, [1], [ind; 62 .- ind])
ind = [4; 16:18; 24:25]
DeltaModel.swap_eigenstates!(H, 68, 73, [2], [ind; 62 .- ind])

colours = [GREY, BLACK, GREEN, BLUE2]
m1 = 64
figa = plot();
for m in m1:m1+11
    for ik in 1:n_cells
        label = ik == 1 && m % 3 == 0 ? "$(m÷3 + skipbands)" : ""
        plot!(figa, φₓ ./ π, H.E[m, ik, :] .+ 8ω .- ω/s*skipbands; label, c=colours[(m-m1)÷3 + 1], xlabel=L"\varphi_x", ylabel="Quasienergy")
    end
end
plot!(xlabel=L"\varphi_t=\varphi_x\ (\pi"*" rad)")
annotate!(figa, 0.5, 240, (L"\nu_2=-1", RED))
annotate!(figa, 0.5, 270, (L"\nu_2=+1", RED))

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
    return plot(figs..., layout=(2, 2), link=:both, plot_title=L"\varphi_t=\varphi_x="*φ_str, xaxis=((0, 2), 0:1:2), yaxis=((0, 2), 0:1:2), legend=false)
end

targetsubbands = [m1+2, m1+11]
DeltaModel.compute_wanniers!(H; targetsubbands)
iφ = 1
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [iφ])
figb = four_wanniers(w, "0", 1:4)

iφ = 11
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [iφ])
figc = four_wanniers(w, L"\pi/6", [3, 4, 2, 1])

plot(figa, figb, figc, layout=grid(1, 3, widths=[0.285, 0.35, 0.35]), plot_title="", dpi=600)
savefig("fig1.png")

##########
########## FIG S2 (Wanniers)
##########

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

set_defaults(width=2*8.6, height=7.5)

DeltaModel.diagonalise!(H, reorder=false)

startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)

DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=true)

iφ = 1
x, _, w = DeltaModel.make_wannierfunctions(H, n_x, Ωt, [iφ])
fig = wanniers12(w)
savefig("figS2.png")

##########
########## FIG S3 (TB Hamiltonian matrix)
##########

set_defaults(width=8.6, height=7.5)

Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)

M = copy(Htb.H[:, :, 1])
M[diagind(M)] .= 0
M = abs.(M)
heatmap(1:size(M, 1), 1:size(M, 1), M; yaxis=:flip, c=CMAP, xticks=1:4:24, yticks=1:4:24, xlabel=L"\ell", ylabel=L"\ell'", cbar_title=L"|J_{\ell'\ell}|", rightmargin=-5mm)
savefig("figS3.pdf")

##########
########## FIG 2
##########

set_defaults(width=2*8.6, height=5)

### (a) Separated spatial spectrum

# spatial part
pumptype = :space
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=false)
startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)
DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=false)
Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)

# spatial part
H_S = DeltaModel.separate_space(Htb; N=n_cells, periodic=true)
S_E = Matrix{Float64}(undef, size(H_S, 1), size(H_S, 3))
for iφ in axes(H_S, 3)
    S_E[:, iφ] = eigvals(Hermitian(H_S[:, :, iφ])) .+ 8ω .- ω/s*skipbands
end
figa = plot();
for r in eachrow(S_E)
    plot!(φₓ ./ π, r, legend=false, xlabel=L"\varphi_x/\pi", ylabel="Quasienergy", c=BLACK, title="Spectrum of "*L"\hat{H}_{{\rm s}}", lw=0.3)
end
plot!(ylims=(245, 270))
annotate!(0.7, 250, (L"\nu_1=+1", RED));
annotate!(0.7, 257, (L"\nu_1=-2", RED))
annotate!(0.7, 267, (L"\nu_1=+1", RED))

### (b) Separated temporal spectrum

pumptype = :time
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=false)
DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=false)
Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)

# T_E = DeltaModel.separate_time_spectrum(Htb, N=2, periodic=true) .+ 4ω
T_E = Matrix{Float64}(undef, 4, size(Htb.H, 3))
for iφ in eachindex(φₓ)
    T_E[:, iφ] = eigvals(Htb.H[1:4, 1:4, iφ])
end
T_E .+= -T_E[1, 1] + (minimum(H.E[:, :, 1]) + 8ω - ω/s*skipbands) - S_E[1, 1]
figb = plot();
for r in eachrow(T_E)
    plot!(φₓ ./ π, r, legend=false, xlabel=L"\varphi_t/\pi", ylabel="Quasienergy", c=BLACK, title="Spectrum of "*L"\hat{H}_{{\rm t}}", lw=0.3)
end
plot!(ylims=(-12.5, 12.5))
annotate!(0.7, -8, (L"\nu_1=-1", RED));
annotate!(0.7, 4, (L"\nu_1=+1", RED))

### (c) HTB' spectrum

sz = size(S_E, 1)*size(T_E, 1)
HTB = Matrix{Float64}(undef, sz, length(φₓ))
for i = axes(S_E, 1)
    for j = axes(T_E, 1)
        HTB[size(T_E, 1)*(i-1)+j, :] .= S_E[i, :] .+ T_E[j, :]
    end
end
figc = plot();
for r in eachrow(HTB)
    plot!(φₓ ./ π, r, legend=false, xlabel=L"\varphi_t=\varphi_x\ (\pi"*" rad)", ylabel="Quasienergy", c=BLACK, title="Spectrum of "*L"\hat{{\cal H}}_{{\rm TB}}^{\prime}", lw=0.3)
end
annotate!(0.7, 240, (L"\nu_2=-1", RED))
annotate!(0.7, 270, (L"\nu_2=+1", RED))

### (d) 4D/8D spectrum

pumptype = :spacetime
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=false)
DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=false)
Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)
DeltaModel.diagonalise!(Htb)
Htb.E .+= 8ω - ω/s*skipbands

# E4D = Matrix{Float64}(undef, 24^2, length(φₓ))
# for i = axes(Htb.E, 1)
#     for j = axes(Htb.E, 1)
#         E4D[24(i-1)+j, :] .= Htb.E[i, :] .+ Htb.E[j, :]
#     end
# end
E8D = Matrix{Float64}(undef, 6, length(φₓ))
E8D[1, :] = 2Htb.E[1, :]
E8D[2, :] = 2Htb.E[4, :]
E8D[3, :] = Htb.E[1, :] .+ Htb.E[5, :]
E8D[4, :] = Htb.E[20, :] .+ Htb.E[24, :]
E8D[5, :] = 2Htb.E[21, :]
E8D[6, :] = 2Htb.E[24, :]

figd = plot();
plot!(φₓ ./ π, E8D[1, :], fillrange=E8D[2, :], fillalpha=0.3, alpha=0, c=BLACK);
plot!(φₓ ./ π, E8D[3, :], fillrange=E8D[4, :], fillalpha=0.3, alpha=0, c=BLACK);
plot!(φₓ ./ π, E8D[5, :], fillrange=E8D[6, :], xlabel=L"\varphi_x=\varphi_y=\varphi_{t_x}=\varphi_{t_y}\ (\pi"*" rad)", ylabel="Quasienergy", c=BLACK, title="Spectrum of "*L"\hat{{\cal H}}_{8{\rm D}}", fillalpha=0.3, alpha=0, legend=false)
hspan!([maximum(E8D[2, :]), minimum(E8D[3, :])], c=RED, alpha=0.5);
hspan!([maximum(E8D[4, :]), minimum(E8D[5, :])], c=RED, alpha=0.5)
annotate!(0.7, 482, (L"\nu_4=-1", RED))
annotate!(0.7, 542, (L"\nu_4=+1", RED))
gap_lo = minimum(E8D[3, :]) - maximum(E8D[2, :])
gap_hi = minimum(E8D[5, :]) - maximum(E8D[4, :])
width_lo = maximum(E8D[2, :]) - minimum(E8D[1, :])
width_hi = maximum(E8D[6, :]) - minimum(E8D[5, :])
gap_lo/width_lo
gap_hi/width_hi

plot(figa, figb, figc, figd, layout=grid(1, 4), plot_title="")
savefig("fig2.pdf")

# sum of two approximate Hamiltonians (does not apper in the paper)

# sort the levels because `HTB` was construced by summing two spectra, and the levels are actually not sorted
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
plot!(φₓ ./ π, E4D[5, :], fillrange=E4D[6, :], xlabel=L"\varphi_x=\varphi_y=\varphi_{t_x}=\varphi_{t_y}\ (\pi"*" rad)", ylabel="Quasienergy", c=BLACK, title="Spectrum of "*L"\hat{{\cal H}}_{8{\rm D}}", fillalpha=0.3, alpha=0, legend=false)
hspan!([maximum(E4D[2, :]), minimum(E4D[3, :])], c=RED, alpha=0.5);
hspan!([maximum(E4D[4, :]), minimum(E4D[5, :])], c=RED, alpha=0.5)
gap_lo = minimum(E4D[3, :]) - maximum(E4D[2, :])
gap_hi = minimum(E4D[5, :]) - maximum(E4D[4, :])
width_lo = maximum(E4D[2, :]) - minimum(E4D[1, :])
width_hi = maximum(E4D[6, :]) - minimum(E4D[5, :])
gap_lo/width_lo
gap_hi/width_hi

##########
########## FIG S4
##########

set_defaults(width=8.6, height=5)

### (a) Spatial pumping

pumptype = :space
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=false)
startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)
DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=true)
Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)

H_S = DeltaModel.separate_space(Htb; N=n_cells, periodic=true)
figa = plot();
for (targetlevels, c) in zip([1:2, 3:4, 5:6], [BLACK, RED, BLUE2])
    pos = DeltaModel.compute_wanniers(H_S; targetlevels, isperiodic=true, n_s=3)
    for r in eachrow(pos)
        scatter!(r, φₓ ./ π; c, legend=false, markerstrokewidth=0, markersize=2)
    end
end
plot!(xlabel=L"i", ylabel=L"\varphi_x/\pi", widen=false, xlims=(0, 6))

### (b) Temporal pumping

pumptype = :time
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=false)
DeltaModel.compute_wanniers!(H; targetsubbands, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=true)
Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)

figb = plot();
for (targetlevels, c) in zip([1:2, 3:4], [BLACK, BLUE2])
    pos = DeltaModel.compute_wanniers(Htb.H[1:4, 1:4, :]; targetlevels, isperiodic=true, n_s=2)
    for r in eachrow(pos)
        scatter!(r, φₓ ./ π; c, legend=false, markerstrokewidth=0, markersize=2)
    end
end
plot!(xlabel=L"\alpha", ylabel=L"\varphi_t/\pi", widen=false, xlims=(0, 4))

plot(figa, figb, layout=(1, 2))
savefig("figS4.pdf")

### (c) Time-space pumping

pumptype = :spacetime
H = DeltaModel.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
DeltaModel.diagonalise!(H, reorder=false)
DeltaModel.compute_wanniers!(H; targetsubbands=1:12, Ωt=0.35π)
DeltaModel.order_wanniers!(H; optimise=true)
Htb = DeltaModel.TBFloquetHamiltonian(H, periodicity=:both)
DeltaModel.diagonalise!(Htb)

targetlevels = 21:24
DeltaModel.compute_wanniers!(Htb, targetlevels)
pos = map(p -> (p/4) divrem(p, 4), Htb.w.pos)
figc = plot()
for p in pos
    scatter!([p[1]], [p[2]], c=:BLACK)
end
figc