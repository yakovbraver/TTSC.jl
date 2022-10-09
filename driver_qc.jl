using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
theme(:dark, size=(800, 600))

includet("bandsolvers.jl")
import .Bandsolvers

include("SpacetimeHamiltonian.jl")

function 𝐻₀(p, x, params)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2
end

function 𝐻(p, x, params, t)
    p^2 + params[1]*cos(2x)^(2params[2]) + params[3]*cos(x)^2 +
    params[4]*sin(2x)^2*cos(2params[6]*t) + 
    params[5]*cos(2x)^2*cos(params[6]*t + pi/2)
end

function 𝑄ₛ(p::Real, x::Real)
    sin(2x)^2
end

function 𝑄ₗ(p::Real, x::Real)
    cos(2x)^2
end

l = 1
gₗ = -7640
Vₗ = -2
λₛ = 100; λₗ = 40; ω = 410
s = 2
params = [gₗ, l, Vₗ, λₛ, λₗ, ω]
H = SpacetimeHamiltonian(𝐻₀, 𝐻, params, s, (1.5, 2), (2, 2.5))

Iₛ, M, coeffs = compute_parameters(H, Function[𝑄ₛ, 𝑄ₗ], [2s, s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
Aₗ = abs(coeffs[2]); χₗ = angle(coeffs[2])

########## Periodic case

φₜ = range(0, 2π, length=61)
n_cells = s
gₗ = -2λₛ*Aₛ
Vₗ = 2λₗ*Aₗ

h = Bandsolvers.UnperturbedHamiltonian(n_cells; M, gₗ, Vₗ, phases=-φₜ/2, maxband=2, isperiodic=true)
Bandsolvers.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(φₜ, r, label=false)
end
plot!(xlabel=L"\phi_t", ylabel="Energy")

# Wannier centres
pyplot()
Bandsolvers.compute_wanniers!(h, targetband=1)
fig = plot();
for (i, ϕ) in enumerate(φₜ)
    scatter!(h.w.pos_lo[i], fill(ϕ, n_cells); marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(h.w.pos_hi[i], fill(ϕ, n_cells); marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\phi_t", cbtitle="Energy")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
w_lo, w_hi = Bandsolvers.make_wannierfunctions(h, x, 1:length(φₜ))
p = Progress(length(φₜ), 1)
@gif for (i, ϕ) in enumerate(φₜ)
    U = @. -λₛ*Aₛ*cos(4x) + λₗ*Aₗ*cos(2x - ϕ) + H.𝐸(Iₛ) - ω/s*Iₛ
    plot(x, U, label=false, ylims=(-5610, -5575))
    scatter!(h.w.pos_lo[i], h.w.E_lo[i]; marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-5610, -5575))
    scatter!(h.w.pos_hi[i], h.w.E_hi[i]; marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:n_cells
        plot!(x, 4abs2.(w_lo[i][j]) .+ h.w.E_lo[i][j], label=false)
        plot!(x, 4abs2.(w_hi[i][j]) .+ h.w.E_hi[i][j], label=false)
    end
    next!(p)
end

########## Non-periodic case

# phases = range(0, π, length=61)
# n_cells = 3
# gₗ = -20; Vₗ = -30
# h = Bandsolvers.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, phases, maxband=2, isperiodic=false)
# Bandsolvers.diagonalise!(h)

# # Energy spectrum
# fig = plot();
# for r in eachrow(h.E)
#     plot!(phases, r, label=false)
# end
# plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))", ylims=(-Inf, 0))
# savefig("nakajima-spectrum.pdf")

# # Wavefunctions
# iϕ = 46; ϕ_str = L"\phi = 3\pi/4"

# x = range(0, n_cells*π, length=100n_cells)
# U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + phases[iϕ])^2
# fig = plot(x ./ π, U, label=false, c=:white, lw=1)

# i = 3 # state number
# ψ = 4abs2.(Bandsolvers.make_eigenfunctions(h, x, [iϕ], [i])) .+ h.E[i, iϕ]
# hline!([h.E[i, iϕ]], c=:white, ls=:dot, lw=0.5, label=false)
# plot!(x ./ π, ψ[:, 1, 1], label=false, title=ϕ_str, xlabel="z", ylabel="Energy")
# savefig("wf-phi-3pi4.pdf")

# # Wannier centres
# n_cells = 3
# gₗ = -20; Vₗ = -30
# phases = [range(0, 0.005, length=10); range(0.006, 3.11, length=40); range(3.14, pi, length=10)]
# h = Bandsolvers.UnperturbedHamiltonian(n_cells; gₗ, Vₗ, phases, maxband=2, isperiodic=false)
# Bandsolvers.diagonalise!(h)
# Bandsolvers.compute_wanniers!(h; targetband=1)

# fig = plot();
# for (i, ϕ) in enumerate(phases)
#     scatter!(h.w.pos_lo[i], fill(ϕ, length(h.w.pos_lo[i])); marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0)
#     scatter!(h.w.pos_hi[i], fill(ϕ, length(h.w.pos_hi[i])); marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
# end
# plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))"*"; non-periodic")
# savefig("nakajima-wannier.pdf")

# x = range(0, n_cells*π, length=50n_cells)
# w_lo, w_hi = Bandsolvers.make_wannierfunctions(h, x, 1:length(phases))
# p = Progress(length(phases), 1)
# @gif for (i, ϕ) in enumerate(phases)
#     U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
#     plot(x, U, label=false, ylims=(-50, 0))
#     scatter!(h.w.pos_lo[i], h.w.E_lo[i]; marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-41, -22))
#     scatter!(h.w.pos_hi[i], h.w.E_hi[i]; marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
#     for j in eachindex(w_lo[i])
#         plot!(x, 4abs2.(w_lo[i][j]) .+ h.w.E_lo[i][j], label=false)
#     end
#     for j in eachindex(w_hi[i])
#         plot!(x, 4abs2.(w_hi[i][j]) .+ h.w.E_hi[i][j], label=false)
#     end
#     next!(p)
# end