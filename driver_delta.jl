using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

includet("DeltaModel.jl")
import .DeltaModel

"Produce an animation of potential if `iϕ == nothing`; otherwise, plot the potential at phase number `iϕ`."
function plot_potential(H; U::Real, U_title::Real=U, lift::Real=0, iϕ::Union{Nothing, <:Real}=nothing)
    (;a, φₓ, N) = H
    x = Float64[]
    barriers = Float64[]
    for i = 0:3N-1
        append!(barriers, [i*a/3 - 0.01, i*a/3 + 0.01])
        append!(x, [i*a/3 + 0.005, (i+1)*a/3 - 0.005])
    end
    append!(barriers, [3N*a/3 - 0.01, 3N*a/3 + 0.01])
    𝑈 = [x -> U * DeltaModel.𝑔(x; n, a) for n = 0:2]
    if iϕ === nothing
        @gif for φ in φₓ
            V = zeros(length(x))
            for n in 1:3
                V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
            end
            plot(x, V.+lift, ylims=(-1.1U, 2U).+lift, lw=2, c=:white, label=false, xlabel=L"x", ylabel="Energy",
                title=L"a=%$a, U=%$U_title, \lambda=%$(uh.λ), \varphi=%$(round(φ, digits=3))", titlepos=:left)
            vspan!(barriers, c=:grey, label=false)
        end
    else
        φ = φₓ[iϕ]
        V = zeros(length(x))
        for n in 1:3
            V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
        end
        plot(x, V.+lift, ylims=(-1.1U, 2U).+lift, lw=2, c=:white, label=false, xlabel=L"x", ylabel="Energy",
            title=L"a=%$a, U=%$U_title, \lambda=%$(500), \varphi=%$(round(φ, digits=3))", titlepos=:left)
        vspan!(barriers, c=:grey, label=false)
    end
end

n_cells = 3
a = 2; λ = 500; U = 3
φₓ = range(0, 2π, length=61)
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, φₓ)

plot_potential(h; U, lift=0)
plot_potential(h; U, iϕ=16)
savefig("potential.pdf")

# dispersion

function plot_dispersion(ε::AbstractVector; φ::Real, uh::DeltaModel.UnperturbedHamiltonian)
    cos_kL = [DeltaModel.cos_ka_tm(E; φ, uh) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(ka)", ylims=(-4, 4), ticks=:native, xlims=(0, ε[end]),
         title=L"U=%$U, a=%$a, \lambda=%$λ, \varphi=%$(round(φ, digits=3))", titlepos=:left, label=false)
    hline!([-1, 1], c=:white, label=false)
    vline!(((1:9) .* pi ./ (uh.a/3)).^2, c=2, label=L"(\frac{\pi n}{a/3})^2", lw=0.5) # analytical energy for a single well of width `a/3`
end

ε = range(U, 2000, length=2000)
plot_dispersion(ε; φ=φₓ[1], uh=h)
xlims!(340, 360)
savefig("dispersion.pdf")

import IntervalRootFinding as iroots
using IntervalArithmetic: (..)

f(E) = DeltaModel.cos_ka(E; φ=0, uh=h)
bounds = (60, 2000)
rts = iroots.roots(f, bounds[1]..bounds[2], iroots.Newton)
z = [rts[i].interval.lo for i in eachindex(rts)]
sort!(z)
scatter!(z, zeros(length(z)))

DeltaModel.diagonalise!(h, length(z), bounds)

# spectrum

fig = plot();
for m in axes(h.E, 2)
    for ik in axes(h.E, 1)
        plot!(φₓ, h.E[ik, m, :])
    end
end
title!(L"N=%$n_cells, a=%$a, U=%$U, \lambda=%$(h.λ)")
ylims!(346.1, 356.2)
savefig("spectrum.pdf")

# eigenfunctions

iϕ = 1
n_x = 50
m = 3
x, ψ = DeltaModel.make_eigenfunctions(h, n_x, m, [iϕ])
fig = plot();
for ik in 1:n_cells
    for b in 1:3
        plot!(x, abs2.(ψ[:, ik, b, 1]) .+ h.E[ik, 3(m-1)+b, iϕ])
    end
end
display(fig)

trapz(f) = ( sum(f) + (f[1] + f[end]) / 2 ) * x[2]-x[1]
trapz(ψ[:, 2, 3, 1] .* conj(ψ[:, 2, 1, 1]))

# Wannier centres

targetband = 3
DeltaModel.compute_wanniers!(h, targetband)

fig = plot();
for (i, ϕ) in enumerate(φₓ)
    for n in 1:n_cells
        scatter!(h.w.pos[:, n, i], fill(ϕ, size(h.w.pos, 1)); marker_z=h.w.E[:, n, i], c=:coolwarm, label=false, markerstrokewidth=0)
    end
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi", cbtitle="Energy",
      title=L"N=%$n_cells, U=%$U, a=%$a, \lambda=%$λ")

# Wannier functions

x, _, w = DeltaModel.make_wannierfunctions(h, n_x, 1:length(φₓ))

n_x = 50
lift = minimum(h.w.E) - 1
lims = (lift-2, maximum(h.w.E)+2)
p = Progress(length(φₓ), 1)
@gif for iϕ in eachindex(φₓ)
    fig = plot_potential(h; U=0.5, U_title=U, lift, iϕ)
    for b in 1:3
        scatter!(h.w.pos[:, b, iϕ], h.w.E[:, b, iϕ]; label=false, markerstrokewidth=0, ylims=lims, c=1:n_cells, markersize=5)
        for j in 1:n_cells
            plot!(x, 0.5abs2.(w[:, j, b, iϕ]) .+ h.w.E[j, b, iϕ], label=false, c=j)
        end
    end
    next!(p)
end

# calculate hopping strength for TB

targetband = 3
d, pos, E = DeltaModel.compute_wanniers(h, targetband)
ws = Matrix{ComplexF64}(undef, 3n_x+1, 3)
fig = plot();
for j in 1:3
    ws[:, j] = sum(d[i, j] * ψ[1:3n_x+1, 1, i, 1] for i = 1:3)
    plot!(x[1:3n_x+1], abs2.(ws[:, j]) .+ E[j])
end
scatter!(pos, E)

J = [d[:, i%3+1]' * (d[:, i] .* h.E[1, range(3(targetband-1)+1, length=3), 1]) for i in 1:3]

###### TB Hamiltonian

n_cells = 3
a = 2; U = 3; 
φₓ = range(0, 2π, length=61)
φₓ = [range(0, pi/3-0.1, length=10); range(pi/3-0.05, pi/3+0.05, length=11);
      range(pi/3+0.1, 4pi/3-0.1, length=25); range(4pi/3-0.05, 4pi/3+0.05, length=11);
      range(4pi/3+0.1, 2pi, length=20)]
htb = DeltaModel.TBHamiltonian(n_cells; a, U, J, isperiodic=true, φₓ)
DeltaModel.diagonalise!(htb)

# Energy spectrum

shift = h.E[1, 3targetband, 1] - htb.E[end, 1]
fig2 = plot();
for r in eachrow(htb.E)
    plot!(φₓ, r .+ shift, xlabel=L"\varphi_x", ylabel="Energy", legend=false)
end
title!("non-periodic TB: "*L"J=%$(J[1]), U=%$U")
savefig("non-periodic.pdf")

# Energy spectrum calculated using k-space representation

ka = [0, pi]
E = DeltaModel.diagonalise_kspace(htb, ka)
fig2 = plot();
for ik in eachindex(ka)
    for m = 1:3
        plot!(φₓ, E[3(ik-1)+m, :], label=L"m = %$m, k=%$(round(ka[ik], digits=3))")
    end
end
title!("TB: "*L"J=%$(J[1]), U=%$U")

# Wannier centres

DeltaModel.compute_wanniers!(htb)
fig2 = plot();
for (iϕ, ϕ) in enumerate(φₓ)
    for b in 1:3
        scatter!((htb.w.pos[:, b, iϕ].+a/6).%(a*n_cells), fill(ϕ, size(htb.w.pos, 1)); marker_z=htb.w.E[:, b, iϕ] .+ shift, c=:coolwarm, label=false, markerstrokewidth=0)
    end
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi", cbtitle="Energy", title="TB: "*L"J=%$(J[1]), U=%$U")

plot(fig, fig2, layout=(2,1))
savefig(fig2, "centres-tb-nonperiodic-8.pdf")

# Wannier functions

wanniers = DeltaModel.make_wannierfunctions(htb, 1:length(φₓ))
htb.w.E .+= shift

lift = minimum(htb.w.E) - 1
lims = (lift-2, maximum(htb.w.E)+2)
x = range(start=a/6, step=a/3, length=3n_cells)
p = Progress(length(φₓ), 1)
@gif for iϕ in eachindex(φₓ)
    fig = plot_potential(htb; U=0.5, U_title=U, lift, iϕ)
    for b in 1:3
        scatter!((htb.w.pos[:, b, iϕ].+a/6).%(a*n_cells), htb.w.E[:, b, iϕ]; label=false, markerstrokewidth=0, ylims=lims, markersize=5, c=1:n_cells)
        for j in 1:n_cells
            plot!(x, abs2.(wanniers[:, j, b, iϕ]) .+ htb.w.E[j, b, iϕ], label=false, c=j)
        end
    end
    next!(p)
end