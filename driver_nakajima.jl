using Plots, LaTeXStrings

include("bandsolvers.jl")

# Energy spectrum

phases = range(0, π, length=61)
n_cells = 5
n_min = 1
n_max = 5
gₗ = -20; Vₗ = -30
e, E = compute_floquet_bands_with_boundary(;n=n_cells, n_min, n_max, phases, s=2, gₗ, Vₗ, λₗ=0, λₛ=0, ω=0, pumptype=:space)

fig = plot();
for r in eachrow(e)
    plot!(phases, r, label=false)
end
plot!(xlabel=L"\phi", ylabel="Energy", title=L"(V_S, V_L) = (20, 30)")
ylims!(-Inf, 0)
savefig("nakajima-spectrum.pdf")

# Wavefunctions

"Reconstruct the coordinate space wavefunction 𝜓(𝑥) = ∑ⱼ𝑐ⱼsin(𝑗𝑥/𝑛) / √(𝑛π/2)"
function make_coordinate_state(x::AbstractVector{<:Real}, coeffs::AbstractVector{<:Number}; n)
    ψ = zeros(eltype(coeffs), length(x))
    for (j, c) in enumerate(coeffs)
        @. ψ += c * sin(j/n * x)
    end
    return ψ ./ sqrt(n*π/2)
end

ϕ = 3pi/4; ϕ_str = L"\phi = 3\pi/4"
ee, EE, c, b = compute_floquet_bands_states(;n=n_cells, n_min, n_max, phases=[ϕ], s=2, gₗ, Vₗ, λₗ=0, λₛ=0, ω=0, pumptype=:space)

x = range(0, n_cells*π, length=2001)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + ϕ)^2
i = 5 # state number
u = 4make_coordinate_state(x, c[:, i], n=n_cells) .+ ee[i]
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
pos_lower, pos_higher, ε_lower, ε_higher = compute_wannier_centres(;N=n_cells, n_target, n_min, n_max, phases, s=2, gₗ=-20, Vₗ=-30, λₗ=0, λₛ=0, ω=0)

fig = plot();
for (i, ϕ) in enumerate(phases)
    scatter!(pos_lower[i],  fill(ϕ, length(pos_lower[i]));  marker_z=ε_lower[i],  c=:coolwarm, label=false, markerstrokewidth=0)
    scatter!(pos_higher[i], fill(ϕ, length(pos_higher[i])); marker_z=ε_higher[i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)")
savefig("nakajima-wannier.pdf")

########## Periodic case

phases = range(0, π, length=61)
n_cells = 4
n_max = 4
gₗ = -20; Vₗ = -30
e, pos_lower, pos_higher, ε_lower, ε_higher = compute_wannier_centres_periodic(; N=n_cells, n_max, n_target=1, phases, gₗ, Vₗ)

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
plot!(minorgrid=true, xlabel=L"z", ylabel=L"\phi", cbtitle="Energy", title=L"(V_S, V_L) = (20, 30)")
savefig("nakajima-wannier-periodic.pdf")