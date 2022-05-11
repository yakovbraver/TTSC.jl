using Plots, LaTeXStrings

include("bandsolvers.jl")

# Energy spectrum

phases = range(0, π, length=61)
n_cells = 5
n_min = 1
n_max = 5
gₗ = -10; Vₗ = -15
e, E = compute_floquet_bands_with_boundary(;n=n_cells, n_min, n_max, phases, s=2, gₗ, Vₗ, λₗ=0, λₛ=0, ω=0, pumptype=:space)

fig = plot();
for r in eachrow(e)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (20, 30)")
ylims!(-Inf, 0)
savefig("nakajima-spectrum.pdf")

# Wavefunctions
ϕ = 3pi/4; ϕ_str = L"\phi = 3\pi/4"
ee, EE, c, b = compute_floquet_bands_states(;n=n_cells, n_min, n_max, phases=[ϕ], s=2, gₗ, Vₗ, λₗ=0, λₛ=0, ω=0, pumptype=:space)

x = range(0, n_cells*π, length=100n_cells)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
i = 5 # state number
u = 4make_sine_state(x, c[:, i], n=n_cells) .+ ee[i]
fig = plot(x ./ π, U, label=false, c=:white, lw=1)
hline!([ee[i]], c=:white, ls=:dot, lw=0.5, label=false)
plot!(x ./ π, u, label=false, title=ϕ_str, xlabel="z", ylabel="Energy")
savefig("wf-phi-3pi4.pdf")

# Wannier centres

phases = [range(0, 0.005, length=10); range(0.006, 3.11, length=101); range(3.14, pi, length=10)]  # values of the adiabatic phase in (S32)
n_cells = 4
n_min = 1
n_max = 5
n_target = 1
pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher = compute_wannier_centres(;N=n_cells, n_target, n_min, n_max, phases, s=2, gₗ, Vₗ, λₗ=0, λₛ=0, ω=0)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[i],  fill(ϕ, length(pos_lower[i]));  marker_z=ε_lower[i],  c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(pos_higher[i], fill(ϕ, length(pos_higher[i])); marker_z=ε_higher[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)")
savefig("nakajima-wannier.pdf")

pyplot()
x = range(0, n_cells*π, length=50n_cells)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(-50, 0))
    scatter!(pos_lower[i],  ε_lower[i]; marker_z=ε_lower[i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims=(-41, -22))
    scatter!(pos_higher[i], ε_higher[i]; marker_z=ε_higher[i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:length(pos_lower[i])
        plot!(x, 4wf_lower[i][:, j] .+ ε_lower[i][j], label=false)
    end
    for j in 1:length(pos_higher[i])
        plot!(x, 4wf_higher[i][:, j] .+ ε_higher[i][j], label=false)
    end
end

########## Periodic case

phases = range(0, π, length=61)
n_cells = 4
n_max = 4
gₗ = -20; Vₗ = -30
e, pos_lower, pos_higher, ε_lower, ε_higher, wf_lower, wf_higher = compute_wannier_centres_periodic(; N=n_cells, n_max, n_target=1, phases, gₗ, Vₗ)

# Energy spectrum
fig = plot();
for r in eachrow(e)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (20, 30)")

# Wannier centres
fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[:, i],  fill(ϕ, n_cells); marker_z=ε_lower[:, i],  c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(pos_higher[:, i], fill(ϕ, n_cells); marker_z=ε_higher[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)"*"; periodic")
savefig("nakajima-wannier-periodic.pdf")

x = range(0, n_cells*π, length=50n_cells)
@gif for (i, ϕ) in enumerate(phases)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
    plot(x, U, label=false, ylims=(-50, 2))
    scatter!(pos_lower[:, i],  ε_lower[:, i]; marker_z=ε_lower[:, i],  c=:coolwarm, label=false,  markerstrokewidth=0, clims=(-41, -22))
    scatter!(pos_higher[:, i], ε_higher[:, i]; marker_z=ε_higher[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
    for j in 1:size(pos_lower, 1)
        plot!(x, 4wf_lower[:, j, i] .+ ε_lower[j, i], label=false)
        plot!(x, 4wf_higher[:, j, i] .+ ε_higher[j, i], label=false)
    end
end