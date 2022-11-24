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
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, isperiodic=true, φₓ)

plot_potential(h; U, lift=0)
plot_potential(h; U, iϕ=16)
savefig("potential.pdf")

# dispersion

function plot_dispersion(ε::AbstractVector; φ::Real, uh::DeltaModel.UnperturbedHamiltonian)
    cos_kL = [DeltaModel.cos_ka(E; φ, uh) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(ka)", ylims=(-4, 4),
         title=L"U=%$U, a=%$a, \lambda=%$λ, \varphi=%$(round(φ, digits=3))", titlepos=:left, label=false)
    hline!([-1, 1], c=:white, label=false)
    vline!(((1:9) .* pi ./ (uh.a/3)).^2, c=2, label=L"(\frac{\pi n}{a/3})^2", lw=0.5) # analytica energy for a single well of width `a/3`
end

ε = range(U, 500, length=2000)
plot_dispersion(ε; φ=φₓ[1], uh=h)
xlims!(340, 360)
savefig("dispersion.pdf")

# spectrum

bandbounds = [(346, 348), (349, 353), (354, 356.5)]
fig = plot();
for (j, bounds) in enumerate(bandbounds)
    display(bounds)
    DeltaModel.diagonalise!(h, bounds)
    for i in 1:n_cells
        plot!(φₓ, h.E[i, :], label=L"m = %$j, k = 2\pi/Na\cdot%$(i-1)", xlabel=L"\varphi_x", ylabel="Energy", legend=:topleft)
    end
end
title!(L"N=%$n_cells, a=%$a, U=%$U, \lambda=%$(h.λ)")
savefig("spectrum.pdf")

# eigenfunctions

DeltaModel.diagonalise!(h, (354, 356.5))
iϕ = 5
n_x = 50
x, ψ = DeltaModel.make_eigenfunctions(h, n_x, [iϕ])
fig = plot();
for n in 1:n_cells
    plot!(x, abs2.(ψ[:, n, 1]) .+ h.E[n, iϕ])
end
ylims!(345, 358)

# Wannier centres

DeltaModel.compute_wanniers!(h)
pyplot()
fig = plot();
m = 3
for (i, ϕ) in enumerate(φₓ)
    scatter!(hs[m].w.pos[:, i], fill(ϕ, size(hs[m].w.pos, 1)); marker_z=hs[m].w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi", cbtitle="Energy",
      title=L"N=%$n_cells, U=%$U, a=%$a, \lambda=%$λ")

# Wannier functions

pyplot()
x, _, w = DeltaModel.make_wannierfunctions(h, n_x, 1:length(φₓ))
lift = minimum(h.w.E) - 0.5
lims = (lift-1, lift+4)
p = Progress(length(φₓ), 1)
@gif for iϕ in eachindex(φₓ)
    fig = plot_potential(h; U=0.2, U_title=U, lift, iϕ)
    scatter!(h.w.pos[:, iϕ], h.w.E[:, iϕ]; label=false, markerstrokewidth=0, ylims=lims, c=1:n_cells, markersize=5)
    for j in 1:n_cells
        plot!(x, abs2.(w[:, j, iϕ]) .+ h.w.E[j, iϕ], label=false, c=j)
    end
    next!(p)
end

n_x = 50
bandbounds = [(346, 348), (349, 353), (354, 356.5)]
hs = [DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, isperiodic=true, φₓ) for _ in 1:3]
ws = [ComplexF64[;;;], ComplexF64[;;;], ComplexF64[;;;]]
x = Float64[]
for (i, bounds) in enumerate(bandbounds)
    DeltaModel.diagonalise!(hs[i], bounds)
    DeltaModel.compute_wanniers!(hs[i])
    x, _, ws[i] = DeltaModel.make_wannierfunctions(hs[i], n_x, 1:length(φₓ))
end
lift = minimum(hs[1].w.E) - 1
lims = (lift-2, maximum(hs[3].w.E)+2)
p = Progress(length(φₓ), 1)
@gif for iϕ in eachindex(φₓ)
    fig = plot_potential(h; U=0.5, U_title=U, lift, iϕ)
    for b in 1:3
        scatter!(hs[b].w.pos[:, iϕ], hs[b].w.E[:, iϕ]; label=false, markerstrokewidth=0, ylims=lims, c=1:n_cells, markersize=5)
        for j in 1:n_cells
            plot!(x, 0.5abs2.(ws[b][:, j, iϕ]) .+ hs[b].w.E[j, iϕ], label=false, c=j)
        end
    end
    next!(p)
end

# calculate hopping strength for TB

n_cells = 1
φₓ = [0]
a = 2; λ = 500; U = 3
n_x = 50
ε = range(U, 500, length=2000)
plot_dispersion(ε; φ=φₓ[1], uh=hs[1])
bandbounds = [(346, 348), (349, 353), (354, 356.5)]
hs = [DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, isperiodic=true, φₓ) for _ in 1:3]

x = Float64[]
ψs = [ComplexF64[;;;], ComplexF64[;;;], ComplexF64[;;;]]
fig = plot();
for (i, bounds) in enumerate(bandbounds)
    display(i)
    DeltaModel.diagonalise!(hs[i], bounds)
    x, ψs[i] = DeltaModel.make_eigenfunctions(hs[i], n_x, [1])
    plot!(x, real.(ψs[i][:, 1, 1]) .+ hs[i].E[1, 1])
end
display(fig)

d, pos, E = DeltaModel.compute_wanniers(hs)
ws = Matrix{ComplexF64}(undef, 3n_cells*n_x+1, 3)
fig = plot();
for j in 1:3
    ws[:, j] = sum(d[k, j] * ψs[k][:, 1, 1] for k = 1:3)
    plot!(x, abs2.(ws[:, j]) .+ E[j])
end
scatter!(pos, E)

J = [d[:, i%3+1]' * (d[:, i] .* [hs[1].E[1, 1], hs[2].E[1, 1], hs[3].E[1, 1]]) for i in 1:3]

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

shift = hs[3].E[1, 1] - htb.E[end, 1]
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
    scatter!((htb.w.pos[:, iϕ].+a/6).%(a*n_cells), fill(ϕ, size(htb.w.pos, 1)); marker_z=htb.w.E[:, iϕ] .+ shift, c=:coolwarm, label=false, markerstrokewidth=0)
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
    scatter!((htb.w.pos[:, iϕ].+a/6).%(a*n_cells), htb.w.E[:, iϕ]; label=false, markerstrokewidth=0, ylims=lims, markersize=5, c=1:3n_cells)
    for j in 1:3n_cells
        plot!(x, abs2.(wanniers[:, j, iϕ]) .+ htb.w.E[j, iϕ], label=false, c=j)
    end
    next!(p)
end