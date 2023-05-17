# A driving script for analysing Hamiltonian (2) from https://arxiv.org/abs/2305.07668
import TTSC.DeltaModel as dm
using Plots, LaTeXStrings, ProgressMeter

plotlyjs()
theme(:dark, size=(800, 500))

"Produce an animation of potential if `iφ == nothing`; otherwise, plot the potential at phase number `iφ`."
function plot_potential(H; U::Real, U_title::Real=U, lift::Real=0, iφ::Union{Nothing, <:Real}=nothing)
    (;a, φₓ, N) = H
    x = Float64[]
    barriers = Float64[]
    for i = 0:3N-1
        append!(barriers, [i*a/3 - 0.01, i*a/3 + 0.01])
        append!(x, [i*a/3 + 0.005, (i+1)*a/3 - 0.005])
    end
    append!(barriers, [3N*a/3 - 0.01, 3N*a/3 + 0.01])
    𝑈 = [x -> U * DeltaModel.𝑔(x; n, a) for n = 0:2]
    if iφ === nothing
        @gif for φ in φₓ
            V = zeros(length(x))
            for n in 1:3
                V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
            end
            plot(x, V.+lift, ylims=(-1.1U, 2U).+lift, lw=2, c=:white, label=false, xlabel=L"x", ylabel="Energy",
                title=L"a=%$a, U=%$U_title, \lambda=%$(h.λ), \varphi=%$(round(φ, digits=3))", titlepos=:left)
            vspan!(barriers, c=:grey, label=false)
        end
    else
        φ = φₓ[iφ]
        V = zeros(length(x))
        for n in 1:3
            V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
        end
        plot(x, V.+lift, ylims=(-1.1U, 2U).+lift, lw=2, c=:white, label=false, xlabel=L"x", ylabel="Energy",
            title=L"a=%$a, U=%$U_title, \lambda=%$(λ), \varphi=%$(round(φ, digits=3))", titlepos=:left)
        vspan!(barriers, c=:grey, label=false)
    end
end

n_cells = 3
a = 4; λ = 2000; U = 7
φₓ = range(0, 2π, length=31)
h = dm.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)

plot_potential(h; U, lift=0)
plot_potential(h; U, iφ=16)
savefig("potential.pdf")

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
savefig("dispersion.pdf")

f(E) = dm.cos_ka(E; φ=0, uh=h)
bounds = (300, 5100)
import IntervalRootFinding as iroots
using IntervalArithmetic: (..)
rts = iroots.roots(f, bounds[1]..bounds[2])
z = [rts[i].interval.lo for i in eachindex(rts)]
sort!(z)
scatter!(z, zeros(length(z)))

dm.diagonalise!(h; bounds)

# spectrum

fig = plot();
for m in axes(h.E, 2)
    for ik in axes(h.E, 1)
        plot!(φₓ, h.E[ik, m, :], xlabel=L"\varphi_x")
    end
end
plot!(legend=false)
savefig("spatial-spectrum.pdf")

# bandgaps

"""
Return gaps between the highest and middle subband of each band in energy spectrum `E`. Number of cells has to be even for correct results.
Assume 2D or 3D system depending on `type` (`2` or `3`).
`iφ` is the phase of minimum of the lowest level of the highest subband.
"""
function get_bandgaps(E; type::Integer, iφ) 
    gaps = Vector{Float64}(undef, size(E, 2)÷3)
    for b in eachindex(gaps)
        h = max(E[1, 3b, 1], E[end, 3b, 1])
        m = min(E[1, 3b, iφ], E[end, 3b, iφ])
        l = max(E[1, 3b-1, iφ], E[end, 3b-1, iφ])
        if type == 2
            gaps[b] = (2m - (h+l)) #/ (2h - 2m)
        else # if type == 3
            gaps[b] = (3m - (l+2h)) #/ (3h - 3m)
        end
    end
    return gaps
end
g = get_bandgaps(h.E, type=2, iφ=6)
fig1 = scatter(h.E[1, 3:3:end, 6], g, xlabel="Level energy", ylabel="2D bandgap", label=L"\lambda=%$λ, U=%$U", markersize=5, markerstrokewidth=0)
savefig("bandgaps.pdf")

# eigenfunctions

iφ = 1
n_x = 50
whichband = 7
x, ψ = dm.make_eigenfunctions(h, n_x, whichband, [iφ])
fig = plot();
for ik in 1:n_cells
    for b in 1:3
        plot!(x, abs2.(ψ[:, ik, b, 1]) .+ h.E[ik, 3(whichband-1)+b, iφ])
    end
end
display(fig)

l = a/3
skipbands = 7 # number of spatial bands that have been skipped by the choice if `bounds` above
sb = 7 # a subband of choice
b = (sb-1)÷3 + 1 + skipbands # band number which subband `sb` belongs to
h.κ[1, 1, sb, 1] / (π/l) # this is approximately `b`

trapz(f) = ( sum(f) - (f[1] + f[end]) / 2 ) * x[2]-x[1]
trapz(ψ[:, 2, 3, 1] .* conj(ψ[:, 2, 1, 1]))

# Wannier centres

targetband = 16
dm.compute_wanniers!(h, targetband)

fig = plot();
for (i, φ) in enumerate(φₓ)
    for b in 1:3
        scatter!(h.w.pos[:, b, i], fill(φ, size(h.w.pos, 1)); marker_z=h.w.E[:, b, i], c=:coolwarm, label=false, markerstrokewidth=0)
    end
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x", cbtitle="Energy", title=L"N=%$n_cells, U=%$U, a=%$a, \lambda=%$λ")

# Wannier functions

n_x = 50
x, _, w = dm.make_wannierfunctions(h, n_x, 1:length(φₓ))

lift = minimum(h.w.E) - 1
lims = (lift-2, maximum(h.w.E)+2)
p = Progress(length(φₓ), 1)
@gif for iφ in eachindex(φₓ)
    fig = plot_potential(h; U=0.5, U_title=U, lift, iφ)
    for b in 1:3
        scatter!(h.w.pos[:, b, iφ], h.w.E[:, b, iφ]; label=false, markerstrokewidth=0, ylims=lims, c=1:n_cells, markersize=5)
        for j in 1:n_cells
            plot!(x, 0.5abs2.(w[:, j, b, iφ]) .+ h.w.E[j, b, iφ], label=false, c=j)
        end
    end
    next!(p)
end

###### TB Hamiltonian with implicit phase

function plot_wanniercentres(htb::dm.AbstractTBHamiltonian)
    plot();
    for (iφ, φ) in enumerate(φₓ)
        for b in 1:3
            scatter!((htb.w.pos[:, b, iφ].+a/6).%(a*htb.N), fill(φ, size(htb.w.pos, 1)); marker_z=htb.w.E[:, b, iφ], c=:coolwarm, label=false, markerstrokewidth=0)
        end
    end
    plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x", cbtitle="Energy")
end

function plot_pumping(htb::dm.AbstractTBHamiltonian)
    wanniers = dm.make_wannierfunctions(htb, 1:length(φₓ))
    lift = minimum(htb.w.E) - 1
    lims = (lift-2, maximum(htb.w.E)+2)
    x = range(start=a/6, step=a/3, length=3htb.N)
    p = Progress(length(φₓ), 1)
    @gif for iφ in eachindex(φₓ)
        fig = plot_potential(htb; U=0.5, U_title=U, lift, iφ)
        for b in 1:3
            scatter!((htb.w.pos[:, b, iφ].+a/6).%(a*htb.N), htb.w.E[:, b, iφ]; label=false, markerstrokewidth=0, ylims=lims, markersize=5, c=1:htb.N)
            for j in 1:htb.N
                plot!(x, abs2.(wanniers[:, j, b, iφ]) .+ htb.w.E[j, b, iφ], label=false, c=j)
            end
        end
        next!(p)
    end
end

# Construct Wanniers at a single phase (any phase is OK)
iφ₀ = 1
targetband = 16
d, pos, E = dm.compute_wanniers(h; targetband, iφ₀)
n_x = 50
x, ψ = dm.make_eigenfunctions(h, n_x, targetband, [iφ₀])
ws = Matrix{ComplexF64}(undef, length(x), 3n_cells)
fig = plot();
for j in 1:3n_cells
    ws[:, j] = sum(d[i, j] * ψ[range((i-1)length(x) + 1, length=length(x))] for i = 1:3n_cells)
    plot!(x, abs2.(ws[:, j]) .+ E[j], c=j)
end
scatter!(pos, E, c=1:3n_cells, markerstrokewidth=0, xlabel=L"x", ylabel="Energy", legend=false)
savefig("tb-wanniers.pdf")

# Construct and diagonalise the TB Hamiltonian
htb = dm.TBHamiltonian(h; d, isperiodic=true, targetband)
dm.diagonalise!(htb)

# Energy spectrum
fig = plot();
for r in eachrow(htb.E)
    plot!(φₓ, r, xlabel=L"\varphi_x", ylabel="Energy", legend=false)
end
plot!(title="TB")
savefig("tb-spectrum.pdf")

plot(φₓ, real(htb.H[1, 1, :]))
plot!(φₓ, real(htb.H[2, 2, :]))
plot!(φₓ, real(htb.H[3, 3, :]))

# Wannier centres
dm.compute_wanniers!(htb)
plot_wanniercentres(htb)

# Pumping animation
plot_pumping(htb)

###### TB Hamiltonian with explicit phase

# calculate hopping strengths between sites 1-2, 2-3, 3-4 at a single phase
iφ = 1
J = [d[:, i+1]' * (d[:, i] .* h.E[range((iφ-1)size(h.E, 2)size(h.E, 1) + 3(targetband-1)size(h.E, 1) + 1, length=3n_cells)]) for i in 1:3]
# construct TB Hamiltonian using only the modulus of 𝐽₁₂
htb = dm.SimpleTBHamiltonian(8; a, U, J=fill(abs(J[1]), 3), isperiodic=false, φₓ)
dm.diagonalise!(htb)

# Energy spectrum
shift = h.E[1, 3targetband, 1] - htb.E[end, 1]
htb.w.E .+= shift # shift to actual energies which are not accounted for in the model
htb.E .+= shift
fig2 = plot();
for r in eachrow(htb.E)
    plot!(φₓ, r, xlabel=L"\varphi_x", ylabel="Energy", legend=false)
end
plot!(title="Explicit TB")
savefig("explicit-tb-spectrum.pdf")

# Energy spectrum calculated using k-space representation

ka = [0, pi]
E = dm.diagonalise_kspace(htb, ka)
fig2 = plot();
for ik in eachindex(ka)
    for m = 1:3
        plot!(φₓ, E[3(ik-1)+m, :], label=L"m = %$m, k=%$(round(ka[ik], digits=3))")
    end
end
title!("TB: "*L"J=%$(abs(J[1])), U=%$U")

# Wannier centres
dm.compute_wanniers!(htb)
plot_wanniercentres(htb)

# Pumping animation
plot_pumping(htb)