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

phases = range(0, π, length=61)
n_cells = 2

h = Bandsolvers.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, phases, maxband=30, isperiodic=true)
Bandsolvers.diagonalise!(h)

# energy of the unperturbed Hamiltonian spectrum
fig = plot();
plot!(range(0, π, length=200), x -> 𝐻₀(0, x, params), lw=2, c=:white, label=false) # spatial potential
for r in eachrow(h.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi_x", ylabel="Energy")

H = Bandsolvers.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:space, minband=1)
Bandsolvers.diagonalise!(H)

# Floquet quasienergy spectrum
fig = plot();
for r in eachrow(H.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Ordered quasienergy spectrum
E_ordered = Bandsolvers.order_floquet_levels(H)
fig = plot();
for (i, r) in enumerate(eachrow(E_ordered))
    n = i + H.minlevel - 1
    b = (n - 1) ÷ 2n_cells + 1
    plot!(phases, r, label="band $b, level $n", c=b)
end
plot!(xlabel=L"\phi_x", ylabel="Quasienergy")

# Maps of Floquet modes
x = range(0, n_cells*pi, length=50n_cells)
Ωt = range(0, 2π, length=40s)
iϕ = 1
whichstates = 1:4
u = Bandsolvers.make_eigenfunctions(H, x, Ωt, [iϕ], whichstates) .|> abs2
figs = [plot() for _ in eachindex(whichstates)]
for (f, n) in enumerate(whichstates)
    figs[f] = heatmap(x, Ωt, u[:, :, n, iϕ]', xlabel=L"x", ylabel=L"\Omega t", c=:viridis, title="Mode $n")
end
plot(figs...)

# # Wannier centres
# pyplot()
# Bandsolvers.compute_wanniers!(h, targetband=1)
# fig = plot();
# for (i, ϕ) in enumerate(φₜ)
#     scatter!(h.w.pos_lo[i], fill(ϕ, n_cells); marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0)
#     scatter!(h.w.pos_hi[i], fill(ϕ, n_cells); marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
# end
# plot!(minorgrid=true, xlabel=L"x", ylabel=L"\phi_t", cbtitle="Energy")

# # Wannier functions
# x = range(0, n_cells*π, length=50n_cells)
# w_lo, w_hi = Bandsolvers.make_wannierfunctions(h, x, 1:length(φₜ))
# p = Progress(length(φₜ), 1)
# @gif for (i, ϕ) in enumerate(φₜ)
#     U = @. -λₛ*Aₛ*cos(4x) + λₗ*Aₗ*cos(2x - ϕ) + H.𝐸(Iₛ) - ω/s*Iₛ
#     plot(x, U, label=false, ylims=(-5610, -5575))
#     scatter!(h.w.pos_lo[i], h.w.E_lo[i]; marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-5610, -5575))
#     scatter!(h.w.pos_hi[i], h.w.E_hi[i]; marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
#     for j in 1:n_cells
#         plot!(x, 4abs2.(w_lo[i][j]) .+ h.w.E_lo[i][j], label=false)
#         plot!(x, 4abs2.(w_hi[i][j]) .+ h.w.E_hi[i][j], label=false)
#     end
#     next!(p)
# end

# ########## Non-periodic case

# h = Bandsolvers.UnperturbedHamiltonian(n_cells; M, gₗ, Vₗ, phases=-φₜ/2, maxband=2, isperiodic=false)
# Bandsolvers.diagonalise!(h)
# h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

# # Energy spectrum
# fig = plot();
# for r in eachrow(h.E)
#     plot!(phases, r, label=false)
# end
# plot!(xlabel=L"\phi_t", ylabel="Energy")

# # Wannier centres
# Bandsolvers.compute_wanniers!(h; targetband=1)

# fig = plot();
# for (i, ϕ) in enumerate(φₜ)
#     scatter!(h.w.pos_lo[i], fill(ϕ, length(h.w.pos_lo[i])); marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0)
#     scatter!(h.w.pos_hi[i], fill(ϕ, length(h.w.pos_hi[i])); marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
# end
# plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi_t", cbtitle="Energy")

# x = range(0, n_cells*π, length=50n_cells)
# w_lo, w_hi = Bandsolvers.make_wannierfunctions(h, x, 1:length(φₜ))
# p = Progress(length(φₜ), 1)
# @gif for (i, ϕ) in enumerate(φₜ)
#     U = @. -λₛ*Aₛ*cos(4x) + λₗ*Aₗ*cos(2x - ϕ) + H.𝐸(Iₛ) - ω/s*Iₛ
#     plot(x, U, label=false, ylims=(-5610, -5575))
#     scatter!(h.w.pos_lo[i], h.w.E_lo[i]; marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-5610, -5575))
#     scatter!(h.w.pos_hi[i], h.w.E_hi[i]; marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
#     for j in eachindex(w_lo[i])
#         plot!(x, 4abs2.(w_lo[i][j]) .+ h.w.E_lo[i][j], label=false)
#     end
#     for j in eachindex(w_hi[i])
#         plot!(x, 4abs2.(w_hi[i][j]) .+ h.w.E_hi[i][j], label=false)
#     end
#     next!(p)
# end