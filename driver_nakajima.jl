using Plots, LaTeXStrings, ProgressMeter
plotlyjs()
theme(:dark, size=(800, 600))

includet("bandsolvers.jl")
import .Bandsolvers

########## Periodic case

phases = [range(0, pi/4-0.1, length=10); range(pi/4-0.01, pi/4+0.01, length=10);
          range(pi/4+0.1, 3pi/4-0.1, length=20); range(3pi/4-0.01, 3pi/4+0.01, length=10);
          range(3pi/4+0.1, pi, length=10)]
n_cells = 3
gₗ = -20; Vₗ = -30

h = Bandsolvers.UnperturbedHamiltonian(n_cells; gₗ, Vₗ, phases, maxband=2, isperiodic=true)
Bandsolvers.diagonalise!(h)

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))", ylims=(-Inf, 0))

# Wannier centres
pyplot()
Bandsolvers.compute_wanniers!(h, targetband=1)
fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(h.w.pos_lo[i], fill(ϕ, n_cells); marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(h.w.pos_hi[i], fill(ϕ, n_cells); marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))"*"; periodic")
savefig("nakajima-wannier-periodic.pdf")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
ψ_lo, ψ_hi = Bandsolvers.make_wannierfunctions(h, x, 1:length(phases))
p = Progress(length(phases), 1)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(-50, 2))
    scatter!(h.w.pos_lo[i], h.w.E_lo[i]; marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-41, -22))
    scatter!(h.w.pos_hi[i], h.w.E_hi[i]; marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:n_cells
        plot!(x, 4abs2.(ψ_lo[i][j]) .+ h.w.E_lo[i][j], label=false)
        plot!(x, 4abs2.(ψ_hi[i][j]) .+ h.w.E_hi[i][j], label=false)
    end
    next!(p)
end

########## Non-periodic case

phases = range(0, π, length=61)
n_cells = 3
gₗ = -20; Vₗ = -30
h = Bandsolvers.UnperturbedHamiltonian(n_cells; gₗ, Vₗ, phases, maxband=2, isperiodic=false)
Bandsolvers.diagonalise!(h)

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))", ylims=(-Inf, 0))
savefig("nakajima-spectrum.pdf")

# Wavefunctions
iϕ = 46; ϕ_str = L"\phi = 3\pi/4"

x = range(0, n_cells*π, length=100n_cells)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + phases[iϕ])^2
fig = plot(x ./ π, U, label=false, c=:white, lw=1)

i = 3 # state number
ψ = 4abs2.(Bandsolvers.make_eigenfunctions(h, x, [iϕ], [i])) .+ h.E[i, iϕ]
hline!([h.E[i, iϕ]], c=:white, ls=:dot, lw=0.5, label=false)
plot!(x ./ π, ψ[:, 1, 1], label=false, title=ϕ_str, xlabel="z", ylabel="Energy")
savefig("wf-phi-3pi4.pdf")

# Wannier centres
n_cells = 3
gₗ = -20; Vₗ = -30
phases = [range(0, 0.005, length=10); range(0.006, 3.11, length=40); range(3.14, pi, length=10)]
h = Bandsolvers.UnperturbedHamiltonian(n_cells; gₗ, Vₗ, phases, maxband=2, isperiodic=false)
Bandsolvers.diagonalise!(h)
Bandsolvers.compute_wanniers!(h; targetband=1)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(h.w.pos_lo[i], fill(ϕ, length(h.w.pos_lo[i])); marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(h.w.pos_hi[i], fill(ϕ, length(h.w.pos_hi[i])); marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (%$(-gₗ), %$(-Vₗ))"*"; non-periodic")
savefig("nakajima-wannier.pdf")

x = range(0, n_cells*π, length=50n_cells)
w_lo, w_hi = Bandsolvers.make_wannierfunctions(h, x, 1:length(phases))
p = Progress(length(phases), 1)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(-50, 0))
    scatter!(h.w.pos_lo[i], h.w.E_lo[i]; marker_z=h.w.E_lo[i], c=:coolwarm, label=false, markerstrokewidth=0, clims=(-41, -22))
    scatter!(h.w.pos_hi[i], h.w.E_hi[i]; marker_z=h.w.E_hi[i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in eachindex(w_lo[i])
        plot!(x, 4abs2.(w_lo[i][j]) .+ h.w.E_lo[i][j], label=false)
    end
    for j in eachindex(w_hi[i])
        plot!(x, 4abs2.(w_hi[i][j]) .+ h.w.E_hi[i][j], label=false)
    end
    next!(p)
end