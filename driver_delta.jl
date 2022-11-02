using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
pyplot()
theme(:dark, size=(800, 500))

includet("DeltaModel.jl")
import .DeltaModel

function plot_potential(uh::DeltaModel.UnperturbedHamiltonian)
    (;a, U, N) = uh
    x = range(0, N * a, length=100)
    𝑈 = [x -> U * DeltaModel.𝑔(x; n, a) for n = 0:2]
    barriers = [-0.01, 0.01]
    for i = 1:3n_cells
        append!(barriers, [i*a/3 - 0.01, i*a/3 + 0.01])
    end
    @gif for φ in uh.φₓ
        V = zeros(length(x))
        for n in 1:3
            V .+= 𝑈[n].(x) .* cos(φ + 2π*(n-1)/3)
        end
        plot(x, V, ylims=(-1.1U, 2U), lw=2, label=false, xlabel=L"x", ylabel="Energy",
             title=L"U=%$U, a=%$a, \varphi=%$(round(φ, digits=3))", titlepos=:left)
        vspan!(barriers, c=:grey, label=false)
    end
end

n_cells = 5
a = 2; λ = 500; U = 3
φₓ = range(0, 2π, length=61)
h = DeltaModel.UnperturbedHamiltonian(n_cells; a, λ, U, isperiodic=true, φₓ)
plot_potential(h)
savefig("potential.pdf")

function plot_dispersion(ε::AbstractVector; φ::Real, uh::DeltaModel.UnperturbedHamiltonian)
    cos_kL = [DeltaModel.cos_ka(E; φ, uh) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(ka)", ylims=(-4, 4),
         title=L"U=%$U, a=%$a, \lambda=%$λ, \varphi=%$(round(φ, digits=3))", titlepos=:left)
    hline!([-1, 1], c=:white, legend=false)
end

ε = range(U, 1000, length=2000)
plot_dispersion(ε; φ=φₓ[6], uh=h)
savefig("solution.html")

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