using Plots, LaTeXStrings
using SpecialFunctions: gamma
using Combinatorics: factorial, binomial
pyplot()
plotlyjs()
theme(:default, size=(700, 600))

### Calculate the spatial Hamiltonian as a function of the action, and compute the associated derivatives

include("bandsolvers.jl")
include("SpacetimeHamiltonian.jl")

function 𝐻₀(p, x, params)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2
end

function 𝐻(p, x, params, t)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2 +
    params[4]*sin(2x)^2*cos(2params[6]*t) + 
    params[5]*cos(2x)^2*cos(params[6]*t + 3pi/2)
end

function 𝑄ₛ(p::Real, x::Real)
    sin(2x)^2
end

function 𝑄ₗ(p::Real, x::Real)
    cos(2x)^2
end

g = 6000; l = 1;
gₗ = -7640 # -2g*factorial(l) / √π / gamma(l + 0.5)
Vₗ = -2
λₛ = 100; λₗ = 40; ω = 410
s = 2
params = [gₗ, l, Vₗ, λₛ, λₗ, ω]
# H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, (0.8, 1.1), (1.2, 1.8), 0.001)
H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, (1.5, 2), (2, 2.5))

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

function plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, φₜ, Iₛ, s, I_min, I_max)
    ϑ = range(0, 2π, length=50)
    I = range(I_min, I_max, length=50)
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    h₀ = H.𝐸(Iₛ) - ω/s*Iₛ
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] + χₛ) + λₗ*Aₗ*cos(s*ϑ[t] + χₗ - φₜ)
    end
    contour(ϑ, I, E', xlabel=L"\Theta"*", rad", ylabel=L"I", cbartitle="Energy \$H\$ (S17)", color=:viridis, minorgrid=true, minorticks=5, levels=30)
    hline!([Iₛ], label=L"I_s = %$(round(Iₛ, sigdigits=4))", c=:white)
    title!(L"\omega = %$ω, s = %$s, M = %$(round(M, sigdigits=2))"*"\n"*
           L"\lambda_L = %$λₗ, A_L = %$(round(Aₗ, sigdigits=2)), \chi_L = %$(round(χₗ, sigdigits=2)),"*"\n"*
           L"\lambda_S = %$λₛ, A_S = %$(round(Aₛ, sigdigits=2)), \chi_S = %$(round(χₛ, sigdigits=2))")
end

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [2s, s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
Aₗ = abs(coeffs[2]); χₗ = angle(coeffs[2])

I_min = 20; I_max = 30
plot_isoenergies(; ω, M, λₛ=λₛ, Aₛ, χₛ, λₗ=λₗ, Aₗ, χₗ, φₜ=π/2, Iₛ, s, I_min, I_max)
savefig("secular-isoenergies.pdf")

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space
pyplot()
fig = plot();
for i in 5:0.2:10
    I, Θ = compute_IΘ(H, i, n_T=100, χ₀=-0.75) # use χ₀ = 0 and ±0.75
    scatter!(Θ, I, xlabel=L"\theta, rad", ylabel=L"I", markerstrokewidth=0, markeralpha=0.6, label=false, minorgrid=true, minorticks=5)
end
ylims!(I_min, I_max); xlims!(0, 2pi)
title!(L"g_l = %$gₗ, V_L = %$Vₗ, \lambda_S = %$λₛ, \lambda_L = %$λₗ, \omega = %$ω")
savefig(fig, "exact-isoenergies.pdf")

### Calculate secular bands

phases = range(0, 2π, length=61) # values of the adiabatic phase in (S32)
n_bands = 4
bands = compute_qc_bands(; n_bands, phases, s, M, λₗAₗ=λₗ*Aₗ, λₛAₛ=λₛ*Aₛ, χₗ, χₛ) .+ H.𝐸(Iₛ) .- ω/s*Iₛ
levels, states = compute_qc_bands_pbc(; n_levels=4, phases, s, M, λₗAₗ=λₗ*Aₗ, λₛAₛ=λₛ*Aₛ, χₗ, χₛ)

fig2 = plot();
for i in 1:n_bands
    plot!(phases, bands[i, :], fillrange=bands[n_bands+i, :], fillalpha=0.35, label="band $i");
end
for lvl in eachrow(levels)
    plot!(phases, lvl .+ H.𝐸(Iₛ) .- ω/s*Iₛ, label=false);
end
xlabel!(L"\varphi_t"*", rad"); ylabel!("Energy of quantised secular "*L"H")
title!(L"\omega = %$ω, M = %$(round(M, sigdigits=2)), \lambda_L = %$λₗ, A_L = %$(round(Aₗ, sigdigits=2)),"*
       L"\lambda_S = %$λₛ, A_S = %$(round(Aₛ, sigdigits=2))")
savefig("semiclassical-bands.pdf")

### Extract tight-binding parameters

function tb_parameters(E_0_0, E_0_pi, E_pi_0)
    J₀ = E_0_0 / 2
    Δ = √(E_0_pi^2 - 4J₀^2)
    ϵ = E_pi_0 / 2J₀
    return J₀, Δ, ϵ
end

gap = bands[1, 1] - bands[2, 1]
w = bands[1, 1] - gap/2

J₀, Δ = tb_parameters(gap/2, bands[1, end÷2]-w)
E0 = @. sqrt(Δ^2*cos(phases)^2 + 4J₀^2)
title!("Fit patameters: "*L"\Delta = %$(round(Δ, sigdigits=3)), J_0 = %$(round(J₀, sigdigits=3)), w = %$(round(w, sigdigits=3))")
plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\cos^{2}\varphi_t+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)
plot!(phases, -E0 .+ w, c=:white, label=false, lw=0.5)

### Calculate Floquet bands
phases = range(0, π, length=61) # values of the adiabatic phase in (S32)
n_min = 1
n_max = 30
n_bands = n_max-n_min+1
eₖ, Eₖ = compute_floquet_bands(;n_min, n_max, phases, s, l, gₗ, Vₗ, λₗ, λₛ, ω, pumptype=:time)
permute_floquet_bands!(Eₖ, eₖ, n_min, ω, s)
fig1 = plot();
for i in 1:2n_bands
    plot!(phases, Eₖ[i, :], fillrange=Eₖ[2n_bands+i, :], fillalpha=0.3, label="m = $(i+2n_min-2)", legend=:outerright)
end
title!(L"V_L = 2, \lambda_S = 40, \lambda_L = 40, \omega = 410")
savefig("2-40-40-410_periodic.pdf")

for (i, m) in enumerate([8, 6, 4, 2, 3, 1, 7, 5, 10, 9])
    plot!(phases, Eₖ[m, :], fillrange=Eₖ[2n_bands+m, :], fillalpha=0.3, label="m = $(i+2n_min-2)", legend=:outerright)
end
for i in 1:2n_bands
    plot!(phases, Eₖ[i, :], fillrange=Eₖ[2n_bands+i, :], fillalpha=0.3, label="")
end

title!("")
ylims!((-5715, -5693))
fig1 = plot();
for i in 1:2n_bands
    plot!(phases, Eₖ[i, :], c=i, label="m = $(i+2n_min-2), k = 0", legend=:outerright)
    plot!(phases, Eₖ[2n_bands+i, :], c=i, ls=:dash, label="m = $(i+2n_min-2), k = 1", legend=:outerright)
end
for i in 1:2n_bands
    plot!(phases, Eₖ[i, :], c=i, label="")
    plot!(phases, Eₖ[2n_bands+i, :], c=i, ls=:dash, label="")
end
display(fig1)

title!("Floquet spectrum, space-time pumping, "*L"V_L=10, \omega=410")
xlabel!(L"2\varphi_t"*", rad")
ylabel!("Floquet quasi-energy "*L"\varepsilon_{k,m}")
xlabel!(L"2\varphi_t=\varphi_x"*", rad")
title!("2D pumping. "*L"\ell = %$l, g = %$g, V_L = %$Vₗ, \lambda_S = %$λₛ, \lambda_L = %$λₗ, \omega = %$ω")
savefig("pumping-spacetime-omega410-V10.pdf")
ylims!(-2760, -2710)
savefig("04-03-Fig12.html")

b = 2
# spatial fit
gap = Eₖ[b+n_bands, 1] - Eₖ[3, 1] |> abs
w = Eₖ[b+n_bands, 1] - gap/2 |> abs
J₀, Δ, ϵ = tb_parameters(Eₖ[b+n_bands, 1]-w, Eₖ[b+n_bands, end÷4]-w, Eₖ[b, 1]-w)
E0 = @. sqrt(Δ^2*sin(phases)^2 + 4J₀^2)
Ek = @. sqrt( Δ^2*sin(phases)^2 + (2J₀*ϵ*cos(phases))^2 )
# plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\sin^{2}2\varphi_x+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)
plot!(phases, E0 .+ w, c=:white, label="BZ centre", legend=:bottomright, lw=0.5)
plot!(phases, Ek .+ w, c=:white, ls=:dash, label="BZ boundary", legend=:bottomright, lw=0.5)
# temporal fit
gap = Eₖ[b, 1] - Eₖ[b+5, 1] |> abs
w = Eₖ[b, 1] - gap/2 |> abs
J₀, Δ = tb_parameters(gap/2, Eₖ[b, end÷4]-w)
E0 = @. sqrt(Δ^2*cos(2phases)^2 + 4J₀^2)
plot!(phases, E0 .+ w, c=:white, label=L"\pm\sqrt{\Delta^{2}\cos^{2}\varphi_t+4J_{0}^{2}}+w", legend=:bottomright, lw=0.5)

plot!(phases, -E0 .+ w, c=:white, label=false, lw=0.5)
plot!(phases, -Ek .+ w, c=:white, ls=:dash, label=false, lw=0.5)
title!("TB fit: "*L"\Delta = %$(round(Δ, sigdigits=3)), J_0 = %$(round(J₀, sigdigits=3)), \epsilon = %$(round(ϵ, sigdigits=3)), w = %$(round(w, sigdigits=6))")
savefig("pumping-space.pdf")

# plot calculated energy spectrum of ℎₖ

fig2 = plot();
plot!(range(0, π, length=200), x -> 𝐻₀(0, x, params), lw=2, c=:white, label=false) # Spatial potential
for i in 1:2n_bands
    plot!(phases, eₖ[i, :], fillrange=eₖ[2n_bands+i, :], fillalpha=0.3)
end
for i in 1:2n_bands
    plot!(phases, eₖ[i, :], label="")
end
title!("Energy spectrum of "*L"h_k"*" (S21), space pumping, "*L"V_L=15")
ylims!(-210, -180)
xlabel!("phi_x")
xlabel!(L"\varphi_x"*", rad"); ylabel!("Eigenenergy "*L"\epsilon_{k,m}"*" of "*L"h_k"*" (S21)")
savefig("h_k-spectrum-V15-zoom+.pdf")

# compare classical vs quantum energies of ℎₖ
scatter(n_min+0.5:n_max+0.5, I -> H.𝐸(I-1), xlabel=L"I=n+1/2", ylabel=L"E", label="classical", legend=:topleft) # not sure why have to use `I-1` instead `I`
scatter!(n_min+0.5:n_max+0.5, eₖ[1:2:2n_bands, 1], label="quantum")

b = 3
shift = Eₖ[1, 1] - bands[1, 1]
plot!(phases./2, bands[1, :].+shift, fillrange=bands[2+1, :].+shift, fillalpha=0.3, label="semiclassical bands 1 and 2", c=:white)
plot!(phases./2, bands[2, :].+shift, fillrange=bands[2+2, :].+shift, fillalpha=0.3, label=false, c=:white)
title!("Temporal pumping, comparison with semiclassical result")
savefig("quantum-vs-semiclassical.pdf")
ylims!(-2740, -2700)
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
title!("4D spacetime bands")
ylims!(-10830, -10815)
savefig("4D-bands.pdf")

### Quasiclassical bands with open boundary conditions

phases = range(0, 2π, length=61)
n_levels = 5
bands, states = compute_qc_bands_obc(; n_levels, phases, M, λₗAₗ=10λₗ*Aₗ, λₛAₛ=2λₛ*Aₛ, χₗ, χₛ, s)

fig = plot();
for i in 1:n_levels
    plot!(phases, bands[i, :], label="")
end
title!("Eigenenergy spectrum of "*L"H"*" (S32) with $n_cells cells and open BC")
xlabel!(L"\varphi_t"*", rad"); ylabel!("Eigenenergy of "*L"H"*" (S32)")
savefig("obc-time-8.pdf")

# plot states

θ = range(0, 2π, length=501)
i_ϕ = 15
U = @. 10λₗ*Aₗ*cos(s*θ - χₗ - phases[i_ϕ]) + λₛ*Aₛ*cos(2s*θ - χₛ) #+ H.𝐸(Iₛ) - ω/s*Iₛ
plot(θ, U, label="potential", c=:white, legend=:outerright)
for i = 1:n_levels
    ψ = 5abs2.(make_sine_state(θ, states[i_ϕ][:, i], n=2)) .+ bands[i, i_ϕ]
    hline!([bands[i, i_ϕ]], c=:white, ls=:dot, label=false); plot!(θ, ψ, label=L"\psi_{%$i}(\theta)")
end
title!("Wavefunctions at "*L"\varphi_t=\pi/2")

@gif for ϕ in phases
    U = @. 10λₗ*Aₗ*cos(s*θ - χₗ - ϕ) + 2λₛ*Aₛ*cos(2s*θ - χₛ)
    plot(θ, U, label=false, ylims=(-150, 150))
end

@gif for i_ϕ in eachindex(phases)
    U = @. (λₗ*Aₗ*cos(s*θ + phases[i_ϕ]) + λₛ*Aₛ*cos(2s*θ))
    plot(θ, U, label="potential", c=:white, legend=:outerright)
    for i = 1:3
        ψ = 5abs2.(make_sine_state(θ, states[i_ϕ][:, i], n=2)) .+ bands[i, i_ϕ]
        hline!([bands[i, i_ϕ]], c=:white, ls=:dot, label=false); plot!(θ, ψ, label=L"\psi_{%$i}(\theta)")
    end
    ylims!(-15, 3)
    xlabel!(L"\theta"*", rad"); ylabel!(L"\psi_n(\theta)")
end

### Floquet bands with open boundary conditions

phases = range(0, π, length=61) # values of the adiabatic phase in (S32)
n_cells = 5
n_min = 24
n_max = 30
e, E = compute_floquet_bands_with_boundary(;n=n_cells, n_min, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, pumptype=:spacetime)
permute_floquet_bands_with_boundary!(E, e; n_cells, n_min, ω, s)

fig = plot();
for r in eachrow(E)
    plot!(phases, r, label=false)
end
title!(L"V_L = 15, \lambda_S = 150, \lambda_L = 55, \omega = 398")
savefig("15-150-55-398.pdf")

### Wannier centre temporal

phases = [range(0, 0.768, length=5); range(0.769, 0.77, length=10); range(0.9, 5.38, length=20); range(5.39, 5.51, length=20); range(5.51, 2pi, length=5)]
phases = range(0, 2pi, length=50);
pos_lower, pos_higher, ε_lower, ε_higher = compute_wannier_centres_qc(; n_levels=10, phases, M, λₗAₗ=10λₗ*Aₗ, λₛAₛ=10λₛ*Aₛ, χₗ, χₛ, s)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[i],  fill(ϕ, length(pos_lower[i]));  marker_z=ε_lower[i],  c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(pos_higher[i], fill(ϕ, length(pos_higher[i])); marker_z=ε_higher[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"\theta", ylabel=L"\phi_t", cbtitle="Quasienergy", title=L"M = -0.082")
savefig("temporal_bad.pdf")

### Wannier centre spatial

phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=200); range(0.9, 2.2, length=20); range(2.3, 2.4, length=200); range(2.4, pi, length=10)]
phases = range(0, pi, length=50)
n_cells = 4
n_max = 15
n_target = 10
pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher = compute_wannier_centres(;N=n_cells, n_target, n_min=1, n_max, phases, s=2, gₗ, Vₗ, λₗ=0, λₛ=0, ω=0)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[i],  fill(ϕ, length(pos_lower[i]));  marker_z=ε_lower[i],  c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(pos_higher[i], fill(ϕ, length(pos_higher[i])); marker_z=ε_higher[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(xlabel=L"\varphi_X", ylabel="Energy", title="Space pumping, band $n_target")
savefig(fig, "10-centres-nonperiodic.pdf")

x = range(0, n_cells*π, length=50n_cells)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(gₗ+Vₗ, 10), xlabel=L"x", ylabel="Energy", title="Space pumping, band $n_target")
    scatter!(pos_lower[i],  ε_lower[i]; marker_z=ε_lower[i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims=(-65, -34))
    scatter!(pos_higher[i], ε_higher[i]; marker_z=ε_higher[i], c=:coolwarm, label=false, markerstrokewidth=0)
    # for j in 1:length(pos_lower[i])
    #     plot!(x, wf_lower[i][:, j] .+ ε_lower[i][j], label=false, c=j)
    # end
    # for j in 1:length(pos_higher[i])
    #     plot!(x, wf_higher[i][:, j] .+ ε_higher[i][j], label=false, c=5+j)
    # end
end

### Periodic
phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=50); range(0.9, 2.2, length=20); range(2.3, 2.4, length=50); range(2.4, pi, length=10)]
phases = range(0, pi, length=61)
n_cells = 4
n_max = 35
n_target = 31
e, pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher = compute_wannier_centres_periodic(; N=n_cells, n_max, n_target, phases, gₗ, Vₗ)
@gif for (i, ϕ) in enumerate(phases)
    U = @. -2000*cos(2x)^2 + 20*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(gₗ+Vₗ, 10))
end
plotlyjs()
fig = plot();
for r in eachrow(e)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\varphi_x", ylabel="Energy", title="band $n_target")
savefig(fig, "$n_target-VL3-spectrum.pdf")

pyplot()
fig = plot();
clims = ( minimum(ε_lower), maximum(ε_higher) )
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[:, i],  fill(ϕ, n_cells); marker_z=ε_lower[:, i],  c=:coolwarm, label=false, markerstrokewidth=0, clims)
    scatter!(pos_higher[:, i], fill(ϕ, n_cells); marker_z=ε_higher[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(xlabel=L"x", ylabel=L"\varphi_x", title="Space pumping, band $n_target")
savefig(fig, "$n_target-centres.pdf")

x = range(0, n_cells*π, length=25n_cells)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + -3*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(400, 500), xlabel=L"x", ylabel="Energy", title="Space pumping, band $n_target")
    scatter!(pos_lower[:, i],  ε_lower[:, i]; marker_z=ε_lower[:, i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims)
    scatter!(pos_higher[:, i], ε_higher[:, i]; marker_z=ε_higher[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:size(pos_lower, 1)
        plot!(x, 4wf_lower[:, j, i] .+ ε_lower[j, i], label=false)
        plot!(x, 4wf_higher[:, j, i] .+ ε_higher[j, i], label=false)
    end
end

# temporal
λₗAₗ=λₗ*Aₗ; λₛAₛ=λₛ*Aₛ
phases = range(0, 2pi, length=61);
θ = range(0, 2π, length=20s)
@gif for (i, ϕ) in enumerate(phases)
    U = @. λₗAₗ*cos(s*θ - χₗ - ϕ) + λₛAₛ*cos(2s*θ - χₛ)
    plot(θ, U, label=false, ylims=(-λₗAₗ-λₛAₛ, λₗAₗ+λₛAₛ))
end

phases = [range(0, 0.768, length=5); range(0.769, 0.77, length=10); range(0.9, 5.38, length=20); range(5.39, 5.51, length=20); range(5.51, 2pi, length=5)]
e, pos_lo, pos_hi, ε_lo, ε_hi, wf_lo, wf_hi = compute_wannier_centres_qc_periodic(; phases, M, λₗAₗ, λₛAₛ, χₗ, χₛ, s)

fig = plot();
for r in eachrow(e)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\varphi_t", ylabel="Energy")

clims = ( minimum(ε_lo), maximum(ε_hi) )
fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lo[:, i],  fill(ϕ, s); marker_z=ε_lo[:, i],  c=:coolwarm, label=false, markerstrokewidth=0, clims)
    scatter!(pos_hi[:, i], fill(ϕ, s); marker_z=ε_hi[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"\theta", ylabel=L"\varphi_t", cbtitle="Energy", title="Time pumping")
savefig(fig, "time-centres.pdf")

@gif for (i, ϕ) in enumerate(phases)
    U = @. λₗAₗ*cos(s*θ - χₗ - ϕ) + λₛAₛ*cos(2s*θ - χₛ)
    plot(θ, U, label=false, ylims=(e[end, length(phases)÷4]-10, λₗAₗ+λₛAₛ), xlabel=L"\theta", ylabel="Energy", title="Time pumping")
    scatter!(pos_lo[:, i],  ε_lo[:, i]; marker_z=ε_lo[:, i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims)
    scatter!(pos_hi[:, i], ε_hi[:, i]; marker_z=ε_hi[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end

ylims = extrema(wf_lo)
i_ϕ = 60
plot(θ, wf_lo[:, 1, i_ϕ]; ylims, label=L"|w_{\alpha=1}(\theta)|^2", palette=palette(:coolwarm, 2))
plot!(θ, wf_lo[:, 2, i_ϕ]; ylims, label=L"|w_{\alpha=2}(\theta)|^2")
plot!(xlabel=L"\theta", ylabel="probability density", title="Quasiclassical Wannier functions; "*L"\varphi_t=2\pi")
savefig("qc-lo-phi=2pi.pdf")

θ = range(0, 2π, length=40s)
@gif for (i, ϕ) in enumerate(phases)
    plot()
    for j in 1:s
        plot!(θ, wf_lo[:, j, i], label=false; ylims)
        # plot!(θ, wf_hi[:, j, i], label=false; ylims)
    end
    plot!(xlabel=L"\theta", title="Lower temporal band, "*L"\varphi_t=%$(round(ϕ, digits=3))")
end

heatmap(θ, phases, wf_hi[:, 1, :]', c=:viridis)
heatmap(θ, phases, wf_hi[:, 2, :]', c=:viridis)
heatmap(θ, phases[1:end-1], (-wf_hi[:, 1, 1:end-1] .+ wf_hi[:, 2, 1:end-1])', c=:coolwarm, xlabel=L"\theta", ylabel=L"\varphi_t", cbar=false)
title!("Quasiclassical Wannier functions; blue: "*L"|w_{\alpha=1}(\theta)|^2"*", red: "*L"|w_{\alpha=2}(\theta)|^2")
savefig("qc-map.pdf")

######## Floquet

phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
phases = range(0, pi, length=61)
n_cells = 1
n_max = 34
n_target = 1
e, E, pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher = compute_floquet_wannier_centres(;N=n_cells, n_target, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, pumptype=:time)

fig = plot();
for r in eachrow(E)
    plot!(2phases, r, label=false)
end
plot!(xlabel=L"\varphi_t=2\varphi_x", ylabel="Quasienergy")
savefig(fig, "timespace-spectrum.pdf")
ylims!(-5716, -5694)

pyplot()
clims = ( minimum(ε_lower), maximum(ε_higher) )
fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[:, i],  fill(2ϕ, n_cells); marker_z=ε_lower[:, i],  c=:coolwarm, label=false, markerstrokewidth=0, clims)
    scatter!(pos_higher[:, i], fill(2ϕ, n_cells); marker_z=ε_higher[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims)
end
plot!(xlabel=L"x", ylabel=L"\varphi_t=2\varphi_x", title=L"\omega t = 0")
savefig(fig, "timespace-centres.pdf")

x = range(0, n_cells*π, length=10n_cells)
@gif for (i, ϕ) in enumerate(phases)
    plot()
    # U = @. gₗ*cos(2x)^2 + -3*cos(x + ϕ)^2
    # plot(x, U, label=false, ylims=(400, 500), xlabel=L"x", ylabel="Energy", title="Space pumping, band $n_target")
    scatter!(pos_lower[:, i],  ε_lower[:, i]; marker_z=ε_lower[:, i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims, ylims=clims.+(-2, 5), xlims=(0, n_cells*π), markersize=7)
    # scatter!(pos_higher[:, i], ε_higher[:, i]; marker_z=ε_higher[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims, xlims=(0, n_cells*π), markersize=7)
    for j in 1:size(pos_lower, 1)
        plot!(x, wf_lower[:, j, i] .+ ε_lower[j, i], label=false, ylims=(-5698, -5695))
        # plot!(x, wf_higher[:, j, i] .+ ε_higher[j, i], label=false, ylims=(-5697, -5694))
    end
    title!("Lower spatial bands, "*L"\omega t = 0, \varphi_t=\varphi_x=%$(round(2ϕ, digits=3))")
end

x = range(0, n_cells*π, length=40n_cells)
ωts = range(0, 2π, length=41) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
clims = extrema(wf_higher)
@gif for (i, ϕ) in enumerate(phases)
    fig1 = heatmap(x, ωts, wf_higher[:, 1, :, i]', xlabel=L"x", ylabel=L"\omega t/s", title=L"|w_{j=1,\beta=1}(x,t)|^2"; c=:viridis, clims, cbar=false)
    fig2 = heatmap(x, ωts, wf_higher[:, 2, :, i]', xlabel=L"x", yformatter=_->"", title=L"|w_{j=1,\beta=2}(x,t)|^2"; c=:viridis, clims)
    plot(fig1, fig2, layout=(1, 2), link=:y, plot_title=L"\varphi_t=%$(round(2ϕ, digits=3))")
end

pyplot()
i_ϕ = 61
fig1 = heatmap(x, ωts, wf_higher[:, 1, :, i_ϕ]', xlabel=L"x", ylabel=L"\omega t/s", title=L"|w_{j=1,\beta=1}(x,t)|^2"; c=:viridis, clims, cbar=false)
fig2 = heatmap(x, ωts, wf_higher[:, 2, :, i_ϕ]', xlabel=L"x", yformatter=_->"", title=L"|w_{j=1,\beta=2}(x,t)|^2"; c=:viridis, clims)
plot(fig1, fig2, layout=(1, 2), link=:y, plot_title=L"\varphi_t=%$(round(2phases[i_ϕ], digits=3))")
savefig("phi=2pi.pdf")