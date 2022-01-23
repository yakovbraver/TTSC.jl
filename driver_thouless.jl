using Plots, LaTeXStrings
using SpecialFunctions: gamma
using Combinatorics: factorial, binomial
pyplot()
plotlyjs()
theme(:dark, size=(700, 600))

### Calculate the spatial Hamiltonian as a function of the action, and compute the associated derivatives

include("SpacetimeHamiltonian.jl")

function 𝐻₀(p, x, params)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2
end

function 𝐻(p, x, params, t)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2 +
    params[4]*sin(2x)^2*cos(2params[6]*t) + 
    params[5]*cos(2x)^2*cos(params[6]*t)
end

function 𝑄ₛ(p::Real, x::Real)
    sin(2x)^2
end

function 𝑄ₗ(p::Real, x::Real)
    cos(2x)^2
end

g = 4200; l = 3;
gₗ = 2g*factorial(l) / √π / gamma(l + 0.5)
Vₗ = 15

λₛ = 360; λₗ = 40; ω = 535.5
s = 2
params = [gₗ, l, Vₗ, λₛ, λₗ, ω]
H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, (0.8, 1.1), (1.2, 1.8), 0.001)

x = range(0, 2π, length=200)
plot!(x, x -> 𝐻₀(0, x, params)) # Spatial potential
surface(x, x, (x, t) -> 𝐻₀(0, x, params) + λₛ*𝑄ₛ(0,x)*cos(2ω*t) + λₗ*𝑄ₗ(0,x)*cos(ω*t), xlabel="x", ylabel="t") # Space-time potential

function plot_actions(H::SpacetimeHamiltonian)
    figs = [plot() for _ in 1:4];
    x = range(0, π, length=200);
    figs[1] = plot(x, H.𝑈, xlabel=L"x", ylabel=L"U(x)=\tilde{g}_\ell\cos^{2\ell}(2x)+V_L\cos^{2}(x)", legend=false);
    title!(L"\ell = %$l, g = %$g, V_L = %$Vₗ");
    I = Dierckx.get_knots(H.𝐸)
    figs[2] = plot(I, H.𝐸(I), xlabel=L"I", ylabel=L"E", legend=false);
    figs[3] = plot(I, H.𝐸′, xlabel=L"I", ylabel=L"dE/dI", legend=false);
    figs[4] = plot(I, H.𝐸″, xlabel=L"I", ylabel=L"d^2E/dI^2", legend=false, ylims=(-30, 30));
    lay = @layout [a{0.5w} grid(3,1)]
    plot(figs..., layout=lay)
    # savefig("H0.pdf")
end

plot_actions(H)
savefig("h_0-parameters.pdf")

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation

function plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, Iₛ, s, I_min)
    ϑ = range(0, 2π, length=50)
    I_max = last(Dierckx.get_knots(H.𝐸))
    I = range(I_min, I_max, length=50)
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    h₀ = H.𝐸(Iₛ) - ω/s*Iₛ
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] + χₛ) + λₗ*Aₗ*cos(s*ϑ[t] + χₗ)
    end
    contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (S17)", color=:viridis, minorgrid=true, minorticks=5)
    hline!([Iₛ], label=L"I_s = %$(round(Iₛ, sigdigits=4))", c=:white)
    title!(L"\omega = %$ω, s = %$s, M = %$(round(M, sigdigits=2))"*"\n"*
           L"\lambda_L = %$λₗ, A_L = %$(round(Aₗ, sigdigits=2)), \chi_L = %$(round(χₗ, sigdigits=2)),"*"\n"*
           L"\lambda_S = %$λₛ, A_S = %$(round(Aₛ, sigdigits=2)), \chi_S = %$(round(χₛ, sigdigits=2))")
end

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [-2s, -s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
ϕₜ = 0.5
eQ = cis(ϕₜ)*coeffs[2]
Aₗ = abs(eQ); χₗ = angle(eQ)

plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, Iₛ, s, I_min=20)
savefig("secular-isoenergies.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

fig = plot();
for i in 30:40
    println(i)
    I, Θ = compute_IΘ(H, i, n_T=150, ϑ₀=0.0)
    scatter!(Θ, I, xlabel=L"\theta, rad", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false, minorgrid=true, minorticks=5)
end
for i in 30:40
    println(i)
    I, Θ = compute_IΘ(H, i, n_T=150, ϑ₀=0.75)
    scatter!(Θ, I, xlabel=L"\theta, rad", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false)
end
ylims!((30, last(Dierckx.get_knots(H.𝐸))))
title!(L"\ell = %$l, g = %$g, V_L = %$Vₗ, \lambda_S = %$λₛ, \lambda_L = %$λₗ, \omega = %$ω")
display(fig)
savefig(fig, "exact-isoenergies.pdf")

### Calculate secular bands

include("bandsolvers.jl")

phases = range(0, π, length=51) # values of the adiabatic phase in (S32)
n_bands = 4
bands = compute_secular_bands(; n_bands, phases, s, M, λₗAₗ=λₗ*Aₗ, λₛAₛ=λₛ*Aₛ) .+ H.𝐸(Iₛ) .- ω/s*Iₛ

fig2 = plot();
for i in 1:n_bands
    plot!(phases, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i");
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Energy of quantised secular "*L"H"*" (S17)")
title!(L"\omega = %$ω, \lambda_L = %$λₗ, \lambda_S = %$λₛ")
savefig("secular-bands.pdf")

### Extract tight-binding parameters

gap = bands[1, 1] - bands[2, 1]
w = bands[1, 1] - gap/2

function tb_parameters(E_0_0, E_0_pi, E_k_pi, k)
    J₀ = E_0_pi / 2
    Δ = √(E_0_0^2 - 4J₀^2)
    # ε = (E_k_pi^2 - Δ^2 - 2J₀^2 * (1+cos(k))) / (2J₀^2 * (1-cos(k))) |> sqrt
    return J₀, Δ#, ε
end

J₀, Δ = tb_parameters(gap/2, bands[1, end÷2]-w, 1.053+w, 1)
E0 = @. sqrt(Δ^2*cos(phases)^2 + 4J₀^2)
title!("Fit patameters: "*L"\Delta = %$(round(Δ, sigdigits=3)), J_0 = %$(round(J₀, sigdigits=3)), w = %$(round(w, sigdigits=3))")
plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\cos^{2}\varphi_t+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)
plot!(phases, -E0 .+ w, c=:white, label=false, lw=0.5)

# k = 1
# Ek = @. Δ^2*sin(phases)^2 + 2J₀^2 * (1+cos(k) + ε^2*sin(phases)^2 * (1-cos(k))) |> sqrt
# plot!(phases, Ek .- w)

### Calculate Floquet bands

phases = range(0, π, length=51) # values of the adiabatic phase in (S32)
n_bands = 45
ee, EE = compute_floquet_bands(;n_bands, phases, s, l, gₗ, Vₗ, λₗ, λₛ, ω, pumptype=:space)

fig1 = plot();
for i in 1:n_bands
    plot!(phases, EE[i, :], fillrange=EE[n_bands+i, :], fillalpha=0.3, label=false)
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Floquet quasi-energy "*L"\varepsilon_{k,m}")
ylims!((EE[end, end÷4], -1500))
ylims!((280, EE[1, end÷4]+10))
ylims!((-1075, -1045))
title!("Pumping in time ("*L"\varphi_x"*" is constant)")
title!("Pumping in space")
title!("Pumping in space ("*L"\varphi_t"*" is constant)")
savefig("pumping-temporal.pdf")

b = 20
# spatial fit
gap = EE[b, 1] - EE[b+2, 1] |> abs
w = EE[b, 1] - gap/2 |> abs
J₀, Δ = tb_parameters(EE[b, end÷4]+w, gap/2, 1.053+w, 1)
E0 = @. sqrt(Δ^2*sin(2phases)^2 + 4J₀^2)
plot!(phases, E0 .- w, c=:white, label=L"\pm\sqrt{\Delta^{2}\sin^{2}2\varphi_x+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)
# temporal fit
gap = EE[b, 1] - EE[b+5, 1] |> abs
w = EE[b, 1] - gap/2 |> abs
J₀, Δ = tb_parameters(gap/2, EE[b, end÷2]+w, 1.053+w, 1)
E0 = @. sqrt(Δ^2*cos(phases)^2 + 4J₀^2)
plot!(phases, E0 .- w, c=:white, label=L"\pm\sqrt{\Delta^{2}\cos^{2}\varphi_t+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)

plot!(phases, -E0 .- w, c=:white, label=false, lw=0.5)
title!("Fit patameters: "*L"\Delta = %$(round(Δ, sigdigits=3)), J_0 = %$(round(J₀, sigdigits=3)), w = %$(round(-w, sigdigits=3))")

fig2 = plot();
for i in 1:2n_bands
    plot!(phases, ee[i, :], fillrange=ee[2n_bands+i, :], fillalpha=0.3, label=false);
end
title!("Energy spectrum of "*L"h_k"*" (S21)")
xlabel!(L"\varphi_x"*", rad"); ylabel!("Eigenenergy "*L"\epsilon_{k,m}"*" of "*L"h_k"*" (S21)")
savefig("h_k-spectrum.pdf")

b = 40
shift = abs(EE[b, 1] - bands[1, 1])
plot!(phases, bands[1, :].-shift, fillrange=bands[4+1, :].-shift, fillalpha=0.3, label="secular bands 1 and 2", c=:white)
plot!(phases, bands[2, :].-shift, fillrange=bands[4+2, :].-shift, fillalpha=0.3, label=false, c=:white)
title!("Pumping in time, comparison with secular result")
savefig("exact-vs-secular.pdf")
findfirst(<(-1050), EE[1:end, 1])
EE[25, 1]

findfirst(>(-1504), EE[1:end, 1])
ee[24*2, 1] - ee[26*2, 1]