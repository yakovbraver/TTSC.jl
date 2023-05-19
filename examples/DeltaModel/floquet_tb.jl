# A driving script for analysing the tight-binding Hamiltonian (4) from https://arxiv.org/abs/2305.07668
import TTSC.DeltaModel as dm
using Plots, LaTeXStrings, ProgressMeter
using LinearAlgebra: diagind

plotlyjs()
theme(:dark, size=(800, 500))
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

# Unperturbed Hamiltonian

n_cells = 2
a = 4; λ = 2000; U = 7
φₓ = range(0, 2pi, length=31)
h = dm.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)

dm.diagonalise!(h; bounds=(300, 10000))

# Floquet Hamiltonian
λₛ = 20; λₗ = 10; ω = 676.8
s = 2
pumptype = :spacetime
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=true)

# Quasienergy spectrum
skipbands = 7 # number of spatial bands that have been skipped by the choice of `bounds` above
fig2 = plot();
for m in axes(h.E, 2)
    for ik in 1:n_cells
        plot!(fig2, φₓ, H.E[m, ik, :] .- ω/s*skipbands, label="$m (b $(H.ν[m])), $ik", c=H.ν[m], xlabel=L"\varphi_x", ylabel="Quasienergy", ticks=:native)
    end
end
plot!(xlabel=L"\varphi_t=\varphi_x", title=L"\omega=%$ω, \lambda_S=%$λₛ, \lambda_L=%$λₗ, \lambda=%$λ, U=%$U")
savefig("$λ-$U-$λₛ-$λₗ-$ω-timespace-spectrum.pdf")

### Floquet TB

# Construct Wanniers

pumptype = :spacetime
H = dm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype)
dm.diagonalise!(H, reorder=false)

startsubband = 1
targetsubbands = range(startsubband, length=12)
n_subbands = length(targetsubbands)

dm.compute_wanniers!(H; targetsubbands, Ωt=0.35pi)
dm.order_wanniers!(H; optimise=true)

theme(:dark, size=(1000, 1000))
iφ = 1
n_x = 50
Ωt = range(0, 2π, length=40s)
x, _, w = dm.make_wannierfunctions(H, n_x, Ωt, [iφ])
figs = [plot() for _ in 1:length(targetsubbands)*n_cells]
for f in eachindex(figs)
    figs[f] = heatmap(x, Ωt, abs2.(w[:, :, f, 1]'), xlabel=L"x", ylabel=L"\Omega t", c=CMAP, title="Wannier $f")
end
plot(figs..., layout=(length(figs)÷4, 4), plot_title=L"\lambda_S=%$λₛ, \lambda_L=%$λₗ, \omega=%$ω, \lambda=%$λ, %$iφ") |> display
savefig("$λ-$U-$λₛ-$λₗ-wanniers-odd-problem.png")

# Construct the TB Floquet Hamiltonian
Htb = dm.TBFloquetHamiltonian(H, periodicity=:both)

# plot the Hamiltonian matrix
M = copy(Htb.H[:, :, 1])
M[diagind(M)] .= 0
M = abs.(M)
heatmap(1:size(M, 1), 1:size(M, 1), M; yaxis=:flip, c=CMAP, xticks=1:size(M, 1), yticks=1:size(M, 1), title=L"\lambda_S=%$λₛ, \lambda_L=%$λₗ, \omega=%$ω, \lambda=%$λ")
fontsize = 7
annotate!([(i, j, Plots.text(round(M[i, j], sigdigits=2), fontsize, RGB(0, 1, 0), :center)) for j in axes(M, 2) for i in axes(M, 1)])
savefig("$λ-$U-$λₛ-$λₗ-HTB.pdf")

plot(φₓ, real(Htb.H[1, 1, :]), label=L"H_{11}")
plot!(φₓ, real(Htb.H[2, 2, :]), label=L"H_{22}", title="Temporal on-site energies", xlabel=L"\varphi_t", ylabel="Quasienergy")
savefig("$λ-$U-$λₛ-$λₗ-temporal-energies.pdf")

plot(φₓ, abs.(Htb.H[1, 2, :]), label=L"H_{12}")
plot!(φₓ, abs.(Htb.H[1, 4, :]), label=L"H_{14}")
plot!(title="Temporal couplings", xlabel=L"\varphi_t", ylabel="Quasienergy")
savefig("$λ-$U-$λₛ-$λₗ-temporal-couplings.pdf")

# Diagonalise TB Hamiltonian
dm.diagonalise!(Htb)
skipbands = 2
fig = plot();
for r in eachrow(Htb.E)
    plot!(φₓ, r .- ω/s*skipbands)
end
plot!(xlabel=L"\varphi_t=\varphi_x", ylabel="Quasienergy", title="Time-space pumping, periodic in both time and space", legend=false)
savefig("$λ-$U-$λₛ-$λₗ-time-spectrum_periodic.pdf")

"Return proper order of states depending on pumping type. (Hardcoded for `length(φₓ) == 31`.)"
function get_order_tb(iφ::Integer; pumptype::Symbol)
    if pumptype == :time
        if iφ < 9
            order = [1, 2, 3, 4]
        elseif iφ < 18
            order = [4, 3, 1, 2]
        else
            order = [2, 1, 4, 3]
        end
    elseif pumptype == :space
        if iφ < 6
            order = [1, 2, 3, 4]
        elseif iφ < 7
            order = [1, 4, 3, 2]
        else
            order = [3, 4, 1, 2]
        end
    else
        if iφ < 6
            order = [1, 2, 3, 4]
        elseif iφ < 16
            order = [3, 4, 1, 2]
        elseif iφ < 20
            order = [1, 2, 4, 3]
        else
            order = [4, 3, 2, 1]
        end
    end
    return order 
end

"Plot TB Wanniers at phase `iφ`. Pass values `iφ` > `length(φₓ)` to go beyound one period."
function plot_tb_wanniers(w, iφ)
    n_t = 4 # number of temporal sites
    figs = [plot() for _ in 1:length(targetlevels)]
    phase = iφ > length(φₓ) ? iφ - length(φₓ) : iφ
    for (f, o) in enumerate(get_order_tb(phase; pumptype))
        if iφ > length(φₓ)
            if pumptype == :space
                f = (f == 1 ? 3 : f == 2 ? 4 : f == 3 ? 1 : 2)
            elseif pumptype == :time
                f = (f == 1 ? 2 : f == 2 ? 1 : f == 3 ? 4 : 3)
            else
                f = (f == 1 ? 4 : f == 4 ? 1 : f == 2 ? 3 : 2)
            end
        end
        figs[f] = heatmap(1:3n_cells, 1:n_t, w[:, :, o, phase]; clims=extrema(w), c=CMAP, xticks=1:3n_cells, yticks=1:n_t,
                          xlabel="spatial site", ylabel="temporal site", title="Wannier $f")
    end
    plot(figs..., plot_title=L"\varphi_x=\varphi_t = %$(round(φₓ[phase]+2pi*(iφ > length(φₓ)), sigdigits=3))", plot_titlelocation=:left)
end

targetlevels = 21:24
dm.compute_wanniers!(Htb, targetlevels)
dm.swap_wanniers!(Htb.w, 2, 3, eachindex(φₓ)) # swap wanniers 2 and 3 so that 1 and 2 (also 3 and 4) occupy the same spatial sites
whichphases = eachindex(φₓ)
w = abs2.(dm.make_wannierfunctions(Htb, whichphases))
plot(Htb.w.pos')

iφ = 1
plot_tb_wanniers(w, iφ)

p = Progress(2length(φₓ), 1)
@gif for iφ in 1:2length(φₓ)
    plot_tb_wanniers(w, iφ)
    next!(p)
end

# Separation
using LinearAlgebra: eigvals, Hermitian

# spatial part
H_S = dm.separate_space(Htb; N=5, periodic=false)
S_E = Matrix{Float64}(undef, size(H_S, 1), size(H_S, 3))
for iφ in axes(H_S, 3)
    S_E[:, iφ] = eigvals(Hermitian(H_S[:, :, iφ])) #.+ 8ω .- ω/s*skipbands
end
figa = plot();
for r in eachrow(S_E)
    plot!(φₓ ./ π, r, legend=false)
end
plot!(xlabel=L"\varphi_x", title="Spectrum of separated spatial lattice")
savefig("$λ-$U-$λₛ-$λₗ-separated-spatial-nonperiodic.pdf")

# temporal part
T_E = Matrix{Float64}(undef, 4, size(Htb.H, 3))
for iφ in eachindex(φₓ)
    T_E[:, iφ] = eigvals(Htb.H[1:4, 1:4, iφ])
end
T_E .-= S_E[1, 1]
fig = plot();
for r in eachrow(T_E)
    plot!(φₓ, r, legend=false)
end
plot!(xlabel=L"\varphi_t", title="Spectrum of separated temporal lattice")
savefig("$λ-$U-$λₛ-$λₗ-separated-temporal_periodic.pdf")