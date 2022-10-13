using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
theme(:dark, size=(800, 600))

includet("bandsolvers.jl")
import .Bandsolvers

########## Periodic case

phases = [range(0, pi/4-0.1, length=10); range(pi/4-0.01, pi/4+0.01, length=10);
          range(pi/4+0.1, 3pi/4-0.1, length=20); range(3pi/4-0.01, 3pi/4+0.01, length=10);
          range(3pi/4+0.1, pi, length=10)]
n_cells = 4
gₗ = -7640; Vₗ = -2  # gₗ = -20; Vₗ = -30 corresponds exactly to the system in Nakajima et al. (https://www.nature.com/articles/nphys3622)

h = Bandsolvers.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, phases, maxband=30, isperiodic=true)
Bandsolvers.diagonalise!(h)

# Energy spectrum
plotlyjs()
fig = plot();
for r in eachrow(h.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))", ylims=(-Inf, 0))

# Wannier centres
pyplot()
Bandsolvers.compute_wanniers!(h, targetband=25)
fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(h.w.pos[:, i], fill(ϕ, 2n_cells); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))"*"; periodic")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
_, w = Bandsolvers.make_wannierfunctions(h, x, 1:length(phases))
lims = (minimum(h.w.E)-0.5, maximum(h.w.E)+2)
p = Progress(length(phases), 1)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=lims)
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims=lims)
    for j in 1:size(w, 2)
        plot!(x, abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false)
    end
    next!(p)
end

########## Non-periodic case

phases = [range(0, 0.005, length=10); range(0.006, 3.11, length=40); range(3.14, pi, length=10)]
n_cells = 4
gₗ = -7640; Vₗ = -2
h = Bandsolvers.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, phases, maxband=30, isperiodic=false)
Bandsolvers.diagonalise!(h)

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))", ylims=(-Inf, 0))

# Wavefunctions
iϕ = 46; ϕ_str = L"\phi = 3\pi/4"

x = range(0, n_cells*π, length=100n_cells)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + phases[iϕ])^2
fig = plot(x ./ π, U, label=false, c=:white, lw=1)

i = 3 # state number
ψ = 4abs2.(Bandsolvers.make_eigenfunctions(h, x, [iϕ], [i])) .+ h.E[i, iϕ]
hline!([h.E[i, iϕ]], c=:white, ls=:dot, lw=0.5, label=false)
plot!(x ./ π, ψ[:, 1, 1], label=false, title=ϕ_str, xlabel="z", ylabel="Energy")

# Wannier centres
Bandsolvers.compute_wanniers!(h, targetband=25)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(h.w.pos[:, i], fill(ϕ, size(h.w.pos, 1)); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))"*"; non-periodic")

x = range(0, n_cells*π, length=50n_cells)
_, w = Bandsolvers.make_wannierfunctions(h, x, 1:length(phases))
lims = (minimum(h.w.E)-0.5, maximum(h.w.E)+2)
p = Progress(length(phases), 1)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=lims)
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims=lims)
    for j in 1:size(w, 2)
        plot!(x, abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false)
    end
    next!(p)
end