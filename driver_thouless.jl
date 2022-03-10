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

g = 6000; l = 1;
gₗ = 2g*factorial(l) / √π / gamma(l + 0.5)
Vₗ = 20

λₛ = 200; λₗ = 75; ω = 391
s = 2
params = [gₗ, l, Vₗ, λₛ, λₗ, ω]
H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, (0.8, 1.1), (1.2, 1.8), 0.001)

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
ϕₜ = 0.0
eQ = cis(ϕₜ)*coeffs[2]
Aₗ = abs(eQ); χₗ = angle(eQ)

I_min = 22
plot_isoenergies(; ω, M, λₛ=λₛ, Aₛ, χₛ, λₗ=λₗ, Aₗ, χₗ, Iₛ, s, I_min)
savefig("secular-isoenergies.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

fig = plot();
for i in 22:0.5:28
    I, Θ = compute_IΘ(H, i, n_T=150, ϑ₀=0.0)
    scatter!(Θ, I, xlabel=L"\theta, rad", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false, minorgrid=true, minorticks=5)
end
for i in 22:0.5:28
    I, Θ = compute_IΘ(H, i, n_T=150, ϑ₀=0.75)
    scatter!(Θ, I, xlabel=L"\theta, rad", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false, minorgrid=true, minorticks=5)
end
ylims!((I_min, last(Dierckx.get_knots(H.𝐸))))
title!(L"\ell = %$l, g = %$g, V_L = %$Vₗ, \lambda_S = %$λₛ, \lambda_L = %$λₗ, \omega = %$ω")
savefig(fig, "exact-isoenergies.pdf")

### Calculate secular bands

include("bandsolvers.jl")

phases = range(0, π, length=51) # values of the adiabatic phase in (S32)
n_bands = 2
bands = compute_qc_bands(; n_bands, phases, s, M, λₗAₗ=λₗ*Aₗ, λₛAₛ=λₛ*Aₛ) .+ H.𝐸(Iₛ) .- ω/s*Iₛ

fig2 = plot();
for i in 1:n_bands
    plot!(phases, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i");
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Energy of quantised secular "*L"H"*" (S17)")
title!(L"M = %$(round(M, sigdigits=2)), \lambda_L = %$λₗ, A_L = %$(round(Aₗ, sigdigits=2)),"*
       L"\lambda_S = %$λₛ, A_S = %$(round(Aₛ, sigdigits=2))")
savefig("semiclassical-bands.pdf")

### Extract tight-binding parameters

function tb_parameters(E_0_0, E_0_pi)
    J₀ = E_0_pi / 2
    Δ = √(E_0_0^2 - 4J₀^2)
    return J₀, Δ
end

gap = bands[1, 1] - bands[2, 1]
w = bands[1, 1] - gap/2

J₀, Δ = tb_parameters(gap/2, bands[1, end÷2]-w)
E0 = @. sqrt(Δ^2*cos(phases)^2 + 4J₀^2)
title!("Fit patameters: "*L"\Delta = %$(round(Δ, sigdigits=3)), J_0 = %$(round(J₀, sigdigits=3)), w = %$(round(w, sigdigits=3))")
plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\cos^{2}\varphi_t+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)
plot!(phases, -E0 .+ w, c=:white, label=false, lw=0.5)

### Calculate Floquet bands

phases = range(0, π, length=51) # values of the adiabatic phase in (S32)
n_min = 50
n_max = 100
n_bands = n_max-n_min+1
eₖ, Eₖ = compute_floquet_bands(;n_min, n_max, phases, s, l, gₗ, Vₗ=Vₗ, λₗ=λₗ, λₛ=λₛ, ω=ω, pumptype=:space)
eₖ, Eₖ = compute_floquet_bands_with_boundary(;n=4, n_min, n_max, phases, s, l, gₗ, Vₗ=Vₗ, λₗ=λₗ, λₛ=λₛ, ω=ω, pumptype=:space)

fig1 = plot();
for i in 1:n_bands
    plot!(phases, Eₖ[i, :], fillrange=Eₖ[n_bands+i, :], fillalpha=0.3, label=false)
end
xlabel!(L"2\varphi_t=\varphi_x"*", rad"); ylabel!("Floquet quasi-energy "*L"\varepsilon_{k,m}")
title!("2D pumping. "*L"\ell = %$l, g = %$g, V_L = %$Vₗ, \lambda_S = %$λₛ, \lambda_L = %$λₗ, \omega = %$ω")
savefig("pumping-both.pdf")
ylims!((4895, 4945))

b = 1
# spatial fit
gap = Eₖ[b, 1] - Eₖ[b+2, 1] |> abs
w = Eₖ[b, 1] - gap/2 |> abs
J₀, Δ = tb_parameters(Eₖ[b, end÷4]-w, gap/2)
E0 = @. sqrt(Δ^2*sin(2phases)^2 + 4J₀^2)
plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\sin^{2}2\varphi_x+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)
# temporal fit
gap = Eₖ[b, 1] - Eₖ[b+5, 1] |> abs
w = Eₖ[b, 1] - gap/2 |> abs
J₀, Δ = tb_parameters(gap/2, Eₖ[b, end÷4]-w)
E0 = @. sqrt(Δ^2*cos(2phases)^2 + 4J₀^2)
plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\cos^{2}\varphi_t+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)

plot!(phases, -E0 .+ w, c=:white, label=false, lw=0.5)
title!("Space pumping. Parameters: "*L"\Delta = %$(round(Δ, sigdigits=3)), J_0 = %$(round(J₀, sigdigits=3)), w = %$(round(-w, sigdigits=3))")
savefig("pumping-space.pdf")

# plot calculated energy spectrum of ℎₖ

fig2 = plot();
x = range(0, π, length=200)
plot!(x, x -> 𝐻₀(0, x, params), lw=2, c=:white, label=false) # Spatial potential
for i in 1:2n_bands
    plot!(phases, eₖ[i, :], fillrange=eₖ[2n_bands+i, :], fillalpha=0.3, label="band $i")
end
for i in 1:n_bands
    plot!(phases, eₖ[i, :], label="")
end
title!("Energy spectrum of "*L"h_k"*" (S21)")
ylims!((7750, 8500))
xlabel!(L"\varphi_x"*", rad"); ylabel!("Eigenenergy "*L"\epsilon_{k,m}"*" of "*L"h_k"*" (S21)")
savefig("h_k-spectrum.pdf")

# compare classical vs quantum energies of ℎₖ
scatter(n_min+0.5:n_max+0.5, I -> H.𝐸(I-1), xlabel=L"I=n+1/2", ylabel=L"E", label="classical", legend=:topleft) # not sure why have to use `I-1` instead `I`
scatter!(n_min+0.5:n_max+0.5, eₖ[1:2:2n_bands, 1], label="quantum")

b = 3
shift = abs(Eₖ[b, 1] - bands[1, 1])
plot!(phases, bands[1, :].-shift, fillrange=bands[2+1, :].-shift, fillalpha=0.3, label="semiclassical bands 1 and 2", c=:white)
plot!(phases, bands[2, :].-shift, fillrange=bands[2+2, :].-shift, fillalpha=0.3, label=false, c=:white)
title!("Pumping in time, comparison with semiclassical result")
savefig("floquet-vs-semiclassical-2.pdf")
findfirst(<(-1390), Eₖ[1:end, 1])
plot(phases, Eₖ[19, :])

### Plot band Minkowski sums

function make_silhouettes(energies, bandnumbers, n_sils)
    simple_bands = Matrix{Float64}(undef, 2n_sils, size(energies, 2))
    n = size(energies, 1) ÷ 2
    for i in 1:n_sils
        simple_bands[i, :] .= max.(energies[bandnumbers[i], :], energies[n+bandnumbers[i], :])
        simple_bands[n_sils+i, :] .= min.(energies[bandnumbers[n_sils+i], :], energies[n+bandnumbers[n_sils+i], :])
    end
    return simple_bands
end

relevant_bands = 1 .+ [0, 2, 1, 5]
relevant_bands = 1 .+ [2, 6, 5, 7]
n_sils = length(relevant_bands) ÷ 2
spacebands = make_silhouettes(Eₖ, relevant_bands, n_sils)

fig1 = plot();
for i in 1:n_sils
    plot!(phases, spacebands[i, :], fillrange=spacebands[n_sils+i, :], fillalpha=0.3, label=false)
end
xlabel!(L"2\varphi_t=\varphi_x"*", rad"); ylabel!("Floquet quasi-energy "*L"\varepsilon_{k,m}")
title!("2D spacetime bands")
savefig("2D-bands.pdf")

function sum_bands(bands)
    n = size(bands, 1)÷2        # number of input bands
    N = round(Int, (n+1)*n÷2)   # number of output bands
    summed = Matrix{Float64}(undef, 2N, size(bands, 2))
    i = 1
    for b1 in 1:n
        for b2 in b1:n
            summed[i, :] = bands[b1, :].+bands[b2, :]
            summed[i+N, :] = bands[b1+n, :].+bands[b2+n, :]
            i += 1 
        end
    end
    summed
end

function plot_summed_bands(bands)
    n = size(bands, 1) ÷ 2
    fig = plot()
    for i in 1:n
        plot!(phases, bands[i, :], fillrange=bands[i+n, :], fillalpha=0.3, label=false)
        hline!([maximum(bands[i, :]), minimum(bands[i+n, :])], c=:white, label=false)
    end
    return fig
end

su = sum_bands(spacebands)
sr = sum_bands(su)
fig = plot_summed_bands(sr)
xlabel!(L"2\varphi_t=\varphi_x"*", rad"); ylabel!("Floquet quasi-energy "*L"\varepsilon_{k,m}")
title!("6D spacetime bands")
ylims!((19720, 19753))
savefig("6D-bands.pdf")

### Boundary conditions

phases = range(0, 2π, length=50)
n_cells = 4
n_bands = 2n_cells + 5
L = -λₗ*Aₗ; S = -λₛ*Aₛ
bands, states = compute_qc_bands_with_boundary(; phases, M=-M, λₗAₗ=L, λₛAₛ=S, n=n_cells)

fig = plot()
for i in 1:n_bands
    plot!(phases, bands[i, :], label="band $i", c=(i > 3 ? i+1 : i))
end
display(fig)
xlabel!("φₜ, rad"); title!("𝑀 increased to 0.5")
xlabel!(L"\varphi_t"*", rad"); ylabel!("Eigenenergy of "*L"H"*" (S32)")
savefig("bands-4-extended.pdf")

# plot states

function make_coordinate_state(x::AbstractVector{<:Real}, coeffs::AbstractVector{<:Number}; n=2)
    ψ = zeros(eltype(coeffs), length(x))
    for (j, c) in enumerate(coeffs)
        @. ψ += c * sin(j/n * x)
    end
    return ψ ./ (n*π/2)
end

x = range(0, n_cells*π, length=101)
i_ϕ = 13
# @gif for i_ϕ = 1:51
    U = @. (L*cos(2x + phases[i_ϕ]) + S*cos(4x))
    # plot(x, U, c=:white, label=false)#, ylims=(-32,32))
    plot(x, U, label="potential", c=:white, legend=:outerright)#, ylims=(-32,32))
# end
for i = 1:n_bands
    ψ = make_coordinate_state(x, states[i_ϕ][:, i], n=n_cells) .+ bands[i, i_ϕ]
    # hline!([bands[i, i_ϕ]], c=:white, ls=:dot, label=false); plot!(x, ψ, label=false, c=(i > 3 ? i+1 : i))
    hline!([bands[i, i_ϕ]], c=:white, ls=:dot, label=false); plot!(x, ψ, label=L"\psi_{%$i}(\theta)", c=(i > 3 ? i+1 : i))
end
xlabel!("θ, rad"); ylabel!("𝜓ₙ")
xlabel!(L"\theta"*", rad"); ylabel!(L"\psi_n(\theta)")
title!("Wavefunctions at φₜ = π/4")
title!("Wavefunctions at "*L"\varphi_t=\pi/2")
savefig("wf-4-extended-pi2.pdf")