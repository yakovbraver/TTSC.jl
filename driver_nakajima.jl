using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
theme(:dark, size=(800, 600))

includet("bandsolvers.jl")

# Energy spectrum

import .Bandsolvers

phases = range(0, π, length=61)
n_cells = 5
n_min = 1
n_max = 5
gₗ = -20; Vₗ = -30
h = Bandsolvers.UnperturbedHamiltonian(n_cells; s=2, gₗ, Vₗ, phases, maxband=2, isperiodic=false)
Bandsolvers.diagonalise!(h)

fig = plot();
for r in eachrow(h.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (20, 30)")
ylims!(-Inf, 0)
savefig("nakajima-spectrum.pdf")

# Wavefunctions
iϕ = 46; ϕ_str = L"\phi = 3\pi/4"

x = range(0, n_cells*π, length=100n_cells)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + phases[iϕ])^2
fig = plot(x ./ π, U, label=false, c=:white, lw=1)

i = 5 # state number
ψ = 4abs2.(Bandsolvers.make_wavefunction(h, [i], [iϕ], x)) .+ h.E[i, iϕ]
hline!([h.E[i, iϕ]], c=:white, ls=:dot, lw=0.5, label=false)
plot!(x ./ π, ψ[:, 1, 1], label=false, title=ϕ_str, xlabel="z", ylabel="Energy")
savefig("wf-phi-3pi4.pdf")

# Wannier centres

phases = [range(0, 0.005, length=10); range(0.006, 3.11, length=101); range(3.14, pi, length=10)]  # values of the adiabatic phase in (S32)
h = Bandsolvers.UnperturbedHamiltonian(n_cells; s=2, gₗ, Vₗ, phases, maxband=2, isperiodic=true)
Bandsolvers.diagonalise!(h)
Bandsolvers.compute_wanniers!(h, targetband=1)
pyplot()
fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(h.w.pos_lo[:, i], fill(ϕ, length(h.w.pos_lo[:, i])); marker_z=h.w.E_lo[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(h.w.pos_hi[:, i], fill(ϕ, length(h.w.pos_hi[:, i])); marker_z=h.w.E_hi[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)")
savefig("nakajima-wannier.pdf")

pyplot()
x = range(0, n_cells*π, length=50n_cells)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(-50, 0))
    scatter!(h.w.pos_lo[:, i],  h.w.E_lo[:, i]; marker_z=h.w.E_lo[:, i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims=(-41, -22))
    scatter!(h.w.pos_hi[:, i], h.w.E_hi[:, i]; marker_z=h.w.E_hi[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:length(h.w.pos_lo[:, i])
        plot!(x, 4wf_lower[i][:, j] .+ h.w.E_lo[:, i][j], label=false)
    end
    for j in 1:length(h.w.pos_hi[:, i])
        plot!(x, 4wf_higher[i][:, j] .+ h.w.E_hi[:, i][j], label=false)
    end
end

########## Periodic case

phases = range(0, π, length=61)
n_cells = 3
n_max = 4
gₗ = -20; Vₗ = -30

h = Bandsolvers.UnperturbedHamiltonian(n_cells; s=2, gₗ, Vₗ, phases, maxband=2, isperiodic=true)
Bandsolvers.diagonalise!(h)
Bandsolvers.compute_wanniers!(h, targetband=1)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(h.w.pos_lo[:, i], fill(ϕ, n_cells); marker_z=h.w.E_lo[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(h.w.pos_hi[:, i], fill(ϕ, n_cells); marker_z=h.w.E_hi[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)"*"; periodic")
savefig("nakajima-wannier-periodic.pdf")

x = range(0, n_cells*π, length=50n_cells)
ψ_lo = 4abs2.(Bandsolvers.make_wavefunction(h, 1:n_cells, 1:length(phases), x, :lo))
ψ_hi = 4abs2.(Bandsolvers.make_wavefunction(h, 1:n_cells, 1:length(phases), x, :hi))

pyplot()
p = Progress(length(phases), 1)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(-50, 2))
    scatter!(h.w.pos_lo[:, i], h.w.E_lo[:, i]; marker_z=h.w.E_lo[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-41, -22))
    scatter!(h.w.pos_hi[:, i], h.w.E_hi[:, i]; marker_z=h.w.E_hi[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:n_cells
        plot!(x, ψ_lo[:, j, i] .+ h.w.E_lo[j, i], label=false)
        plot!(x, ψ_hi[:, j, i] .+ h.w.E_hi[j, i], label=false)
    end
    next!(p)
end