using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

includet("DeltaModel.jl")
import .DeltaModel

"Produce an animation of potential if `iϕ == nothing`; otherwise, plot the potential at phase number `iϕ`."
function plot_potential(uh::DeltaModel.UnperturbedHamiltonian; U::Real, lift::Real=0, iϕ::Union{Nothing, <:Real}=nothing)
    (;a, φₓ, N) = uh
    x = range(0, N * a, length=100)
    𝑈 = [x -> U * DeltaModel.𝑔(x; n, a) for n = 0:2]
    barriers = [-0.01, 0.01]
    for i = 1:3n_cells
        append!(barriers, [i*a/3 - 0.01, i*a/3 + 0.01])
    end
    if iϕ === nothing
        @gif for φ in φₓ
            V = zeros(length(x))
            for n in 1:3
                V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
            end
            plot(x, V.+lift, ylims=(-1.1U, 2U).+lift, lw=2, c=:white, label=false, xlabel=L"x", ylabel="Energy",
                title=L"U=%$U, a=%$a, \varphi=%$(round(φ, digits=3))", titlepos=:left)
            vspan!(barriers, c=:grey, label=false)
        end
    else
        φ = φₓ[iϕ]
        V = zeros(length(x))
        for n in 1:3
            V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
        end
        plot(x, V.+lift, ylims=(-1.1U, 2U).+lift, lw=2, c=:white, label=false, xlabel=L"x", ylabel="Energy",
            title=L"U=%$U, a=%$a, \varphi=%$(round(φ, digits=3))", titlepos=:left)
        vspan!(barriers, c=:grey, label=false)
    end
end

n_cells = 5
a = 2; λ = 500; U = 3
φₓ = range(0, 2π, length=61)
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, isperiodic=true, φₓ)

plot_potential(h; U, lift=0)
savefig("potential.pdf")

# dispersion

function plot_dispersion(ε::AbstractVector; φ::Real, uh::DeltaModel.UnperturbedHamiltonian)
    cos_kL = [DeltaModel.cos_ka(E; φ, uh) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(ka)", ylims=(-4, 4),
         title=L"U=%$U, a=%$a, \lambda=%$λ, \varphi=%$(round(φ, digits=3))", titlepos=:left)
    hline!([-1, 1], c=:white, legend=false)
end

ε = range(U, 1000, length=2000)
plot_dispersion(ε; φ=φₓ[6], uh=h)
savefig("solution.html")

# spectrum

bandbounds = [(346, 348), (349, 353), (354, 356.5)]
fig = plot()
for (j, bounds) in enumerate(bandbounds)
    display(bounds)
    DeltaModel.diagonalise!(h, bounds)
    for i in 1:n_cells
        plot!(φₓ, h.E[i, :], label="m = $j, k = 2π/Na*$(i-1)", xlabel=L"\varphi_x", ylabel="Energy")
    end
end
display(fig)

# eigenfunctions

DeltaModel.diagonalise!(h, (354, 356.5))
iϕ = 5
n_x = 50
x, ψ = DeltaModel.make_eigenfunctions(h, n_x, [iϕ])
fig = plot()
for n in 1:n_cells
    plot!(x, abs2.(ψ[:, n, 1]) .+ h.E[n, iϕ])
end
ylims!(h.E[1, iϕ]-5, h.E[1, iϕ]+5)

# Wannier centres

DeltaModel.compute_wanniers!(h)

fig = plot();
for (i, ϕ) in enumerate(φₓ)
    scatter!(h.w.pos[:, i], fill(ϕ, size(h.w.pos, 1)); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi", cbtitle="Energy")

# Wannier functions

x, _, w = DeltaModel.make_wannierfunctions(h, n_x, 1:length(φₓ))
lift = minimum(h.w.E) - 0.5
lims = (lift-1, lift+4)
p = Progress(length(φₓ), 1)
@gif for iϕ in eachindex(φₓ)
    fig = plot_potential(h; U=0.2, lift, iϕ)
    scatter!(h.w.pos[:, iϕ], h.w.E[:, iϕ]; label=false, markerstrokewidth=0, ylims=lims, c=1:n_cells, markersize=5)
    for j in 1:n_cells
        plot!(x, abs2.(w[:, j, iϕ]) .+ h.w.E[j, iϕ], label=false, c=j)
    end
    next!(p)
end