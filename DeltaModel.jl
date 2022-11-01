using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
# pyplot()
theme(:dark, size=(800, 500))

function 𝑔(x; n, L)
    Int( n/3 <= (x % L)/L < (n+1)/3 )
end

function plot_potential(;L, U, n_cells)
    x = range(0, n_cells * L, length=100)
    𝑈 = [x -> U * 𝑔(x; n, L) for n = 0:2]
    # ϕ = 0
    barriers = [-0.01, 0.01]
    for i = 1:3n_cells
        append!(barriers, [i/3 - 0.01, i/3 + 0.01])
    end
    @gif for ϕ in range(0, 2π, 61)
        V = zeros(length(x))
        for n in 1:3
            V .+= 𝑈[n].(x) .* cos(ϕ + 2π*(n-1)/3)
        end
        plot(x, V, ylims=(-1.1U, 2U), lw=2, label=false, xlabel=L"x", ylabel="Energy",
             title=L"U=%$U, L=%$L, \phi=%$(round(ϕ, digits=3))", titlepos=:left)
        vspan!(barriers, c=:grey, label=false)
    end
end

plot_potential(L=1, U=1, n_cells=2)
savefig("potential.pdf")

function K(ε; ϕ, L, U, λ)
    κ = [sqrt(ε - U*cos(ϕ + 2π*n/3)) for n = 0:2]
    s = sin.(κ*L/3)
    c = cos.(κ*L/3)
    return (-κ[1]^2*s[1] * (s[2]*(λ*s[3]+κ[3]c[3]) + κ[2]s[3]c[2]) +
            κ[1]c[1] * (s[2] * (3λ*κ[3]c[3] - s[3]*(κ[2]^2+κ[3]^2-2λ^2)) + κ[2]c[2]*(3λ*s[3]+2κ[3]*c[3])) +
            s[1] * (s[2]*(λ*s[3]*(-κ[2]^2-κ[3]^2+λ^2) + κ[3]c[3]*(2λ^2-κ[2]^2)) + κ[2]c[2]*(3λ*κ[3]c[3] - s[3]*(κ[3]^2-2λ^2)))
    ) / (2κ[1]κ[2]κ[3])
end

function plot_dispersion(ε; ϕ, L, U, λ)
    cos_kL = [K(E; ϕ, L, U, λ) for E in ε]
    plot(ε, cos_kL, xlabel=L"\varepsilon", ylabel=L"\cos(kL)", ylims=(-4, 4),
         title=L"U=%$U, L=%$L, \lambda=%$λ, \phi=%$(round(ϕ, digits=3))", titlepos=:left)
    hline!([-1, 1], c=:white, legend=false)
end

U = 3
ε = range(U, 1000, length=2000)
plot_dispersion(ε; ϕ=0, L=2, U, λ=500)
savefig("solution.html")