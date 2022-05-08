using Plots, LaTeXStrings

include("bandsolvers.jl")

phases = range(0, π, length=61) # values of the adiabatic phase in (S32)
n_cells = 5
n_min = 1
n_max = 5
e, E = compute_floquet_bands_with_boundary(;n=n_cells, n_min, n_max, phases, s=2, gₗ=-20, Vₗ=-30, λₗ=0, λₛ=0, ω=0, pumptype=:space)

fig = plot();
for r in eachrow(e)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (20, 30)")
ylims!(-Inf, 0)
savefig("nakajima-spectrum.pdf")

### Wannier centres

phases = [range(0, 0.005, length=10); range(0.006, 3.11, length=101); range(3.14, pi, length=10)]  # values of the adiabatic phase in (S32)
n_cells = 4
n_min = 1
n_max = 5
n_target = 1
pos_lower, pos_higher, ε_lower, ε_higher = compute_wannier_centres(;N=n_cells, n_target, n_min, n_max, phases, s=2, gₗ=-20, Vₗ=-30, λₗ=0, λₛ=0, ω=0)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[i], fill(ϕ, length(pos_lower[i])); marker_z=ε_lower[i],     c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(pos_higher[i], fill(ϕ, length(pos_higher[i])); marker_z=ε_higher[i],  c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)")
savefig("nakajima-wannier.pdf")