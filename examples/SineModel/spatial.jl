# A driving script for analysing unperturbed spatial Hamiltonian (2) from https://doi.org/10.1103/PhysRevB.106.144301 (https://arxiv.org/abs/2206.14804)
import TTSC.SineModel as sm
using Plots, LaTeXStrings, ProgressMeter

plotlyjs()
theme(:dark, size=(800, 500))

using LinearAlgebra.BLAS: set_num_threads
set_num_threads(1)

########## Periodic case

φₓ = [range(0, pi/4-0.1, length=10); range(pi/4-0.01, pi/4+0.01, length=10);
      range(pi/4+0.1, 3pi/4-0.1, length=20); range(3pi/4-0.01, 3pi/4+0.01, length=10);
      range(3pi/4+0.1, pi, length=10)]
n_cells = 4
gₗ = -7640; Vₗ = -2
# gₗ = -30; Vₗ = -20 # corresponds exactly to the system in Nakajima et al. (https://www.nature.com/articles/nphys3622)

h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=32, isperiodic=true)
sm.diagonalise!(h)

# Energy spectrum
fig = plot();
for (i, r) in enumerate(eachrow(h.E))
    plot!(φₓ, r, label="$i")
end
plot!(xlabel=L"\varphi_x", ylabel="Energy")

# Eigenfunctions
iφ = 1;
x = range(0, n_cells*π, length=1000n_cells)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + φₓ[iφ])^2
fig = plot(x ./ π, U, label=false, c=:white, lw=1)
for i in 145:150
    ψ = 4abs2.(sm.make_eigenfunctions(h, x, [iφ], [i])) .+ h.E[i, iφ]
    hline!([h.E[i, iφ]], c=:white, ls=:dot, lw=0.5, label=false)
    plot!(x ./ π, ψ[:, 1, 1], xlabel=L"x", ylabel="Energy", c=i, label="$i")
end
display(fig)

# Wannier centres
sm.compute_wanniers!(h, targetband=27, mixsubbands=false)
fig = plot();
for (i, φ) in enumerate(φₓ)
    scatter!(h.w.pos[:, i], fill(φ, 2n_cells); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x", cbtitle="Energy")

# Wannier functions
x = range(0, n_cells*π, length=500n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₓ))
x_U = range(0, n_cells*π, length=50n_cells)
lims = (minimum(h.w.E)-0.5, maximum(h.w.E)+2)
p = Progress(length(φₓ), 1)
@gif for (i, φ) in enumerate(φₓ)
    U = @. gₗ*cos(2x_U)^2 + Vₗ*cos(x_U + φ)^2
    plot(x_U ./ π, U, label=false, ylims=lims, c=:white)
    scatter!(h.w.pos[:, i] ./ π, h.w.E[:, i]; label=false, markerstrokewidth=0, c=1:2n_cells, markersize=5)
    for j in axes(w, 2)
        plot!(x ./ π, abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false, c=j)
    end
    plot!(xlabel=L"x/\pi", ylabel="Energy", title=L"V_{{\rm S}} = %$gₗ, V_{{\rm L}} = %$Vₗ; \varphi_x = %$(round(φ, sigdigits=3))", titlepos=:left, xlims=(0, n_cells))
    next!(p)
end

##### Tight-binding 

# Compute Wanniers by mixing the two subbands of `targetband` together. Wanniers are construced at all phases,
# even though we need them at only one phase when constructing the TB Hamiltonian.
sm.compute_wanniers!(h, targetband=27, mixsubbands=true)

# Plot the Wanniers
x = range(0, n_cells*π, length=50n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₓ))
lims = (minimum(h.w.E)-0.5, maximum(h.w.E)+2)
iφ = 2
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + φₓ[iφ])^2
fig = plot(x, U, label=false, ylims=lims, c=:white);
scatter!(h.w.pos[:, iφ], h.w.E[:, iφ], label=false, markerstrokewidth=0, c=1:2n_cells);
for j in axes(w, 2)
    plot!(x, abs2.(w[:, j, iφ]) .+ h.w.E[j, iφ], label=false, c=j)
end
display(fig)

# Construct the TB Hamiltonian
htb = sm.TBHamiltonian(h)
sm.diagonalise!(htb)

fig = plot();
for r in eachrow(htb.E)
    plot!(φₓ, r)
end
plot!(xlabel=L"\varphi_x", ylabel="Energy")

# Wannier centres
sm.compute_wanniers!(htb)
fig2 = plot();
for (iφ, φ) in enumerate(φₓ)
    scatter!(htb.w.pos[:, iφ], fill(φ, size(htb.w.pos, 1)); marker_z=htb.w.E[:, iφ], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x", cbtitle="Energy")

# Wannier functions
whichphases = 1:length(φₓ)
wanniers = sm.make_wannierfunctions(htb, whichphases)

lims = (minimum(htb.w.E)-1, maximum(htb.w.E)+1)
x_U = range(0, n_cells*π, length=1000n_cells)
x = range(start=0, step=π/2, length=2n_cells)
p = Progress(length(φₓ), 1)
@gif for (i, iφ) in enumerate(whichphases)
    U = @. gₗ*cos(2x_U)^2 + Vₗ*cos(x_U + φₓ[iφ])^2
    plot(x_U ./ π, U, label=false, ylims=lims, c=:white);
    scatter!(htb.w.pos[:, iφ] ./ π, htb.w.E[:, iφ]; label=false, markerstrokewidth=0, ylims=lims, markersize=5, c=1:2n_cells, xlims=(0, n_cells))
    for j in axes(wanniers, 2)
        plot!(x ./ π, abs2.(wanniers[:, j, i]) .+ htb.w.E[j, iφ], label=false, c=j)
    end
    plot!(xlabel=L"x/\pi", ylabel="Energy", title=L"V_{{\rm S}} = %$gₗ, V_{{\rm L}} = %$Vₗ; \varphi_x = %$(round(φₓ[iφ], sigdigits=3))", titlepos=:left)
    next!(p)
end

########## Non-periodic case

# φₓ = [range(0, 0.005, length=10); range(0.006, 3.11, length=40); range(3.14, pi, length=10)] # good for the system of Nakajima et al.
n_cells = 4
gₗ = -7640; Vₗ = -2
h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=32, isperiodic=false)
sm.diagonalise!(h)

# Energy spectrum
fig = plot();
for r in eachrow(h.E)
    plot!(φₓ, r, label=false)
end
plot!(xlabel=L"\varphi_x", ylabel="Energy", ylims=(-Inf, 0))

# Wavefunctions
iφ = 46; φ_str = L"\varphi_x = 3\pi/4"

x = range(0, n_cells*π, length=100n_cells)
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + φₓ[iφ])^2
fig = plot(x ./ π, U, label=false, c=:white, lw=1)

i = 3 # state number
ψ = 4abs2.(sm.make_eigenfunctions(h, x, [iφ], [i])) .+ h.E[i, iφ]
hline!([h.E[i, iφ]], c=:white, ls=:dot, lw=0.5, label=false)
plot!(x ./ π, ψ[:, 1, 1], label=false, title=φ_str, xlabel="z", ylabel="Energy")

# Wannier centres
sm.compute_wanniers!(h, targetband=25, mixsubbands=false)

fig = plot();
for (i, φ) in enumerate(φₓ)
    scatter!(h.w.pos[:, i], fill(φ, size(h.w.pos, 1)); marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0)
end
plot!(minorgrid=true, xlabel=L"x", ylabel=L"\varphi_x", cbtitle="Energy")

# Wannier functions
x = range(0, n_cells*π, length=50n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₓ))
lims = (minimum(h.w.E)-0.5, maximum(h.w.E)+2)
p = Progress(length(φₓ), 1)
@gif for (i, φ) in enumerate(φₓ)
    U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + φ)^2
    plot(x, U, label=false, ylims=lims)
    scatter!(h.w.pos[:, i], h.w.E[:, i]; marker_z=h.w.E[:, i], c=:coolwarm, label=false, markerstrokewidth=0, clims=lims)
    for j in axes(w, 2)
        plot!(x, abs2.(w[:, j, i]) .+ h.w.E[j, i], label=false)
    end
    next!(p)
end

##### Wannier functions obtained by mixing all subbands

sm.compute_wanniers!(h, targetband=25, mixsubbands=true)

x = range(0, n_cells*π, length=50n_cells)
_, w = sm.make_wannierfunctions(h, x, 1:length(φₓ))
lims = (minimum(h.w.E)-0.5, maximum(h.w.E)+2)
iφ = 1
U = @. gₗ*cos(2x)^2 + Vₗ*cos(x + φₓ[iφ])^2
plot(x, U, label=false, ylims=lims);
scatter!(h.w.pos[:, iφ], h.w.E[:, iφ]; marker_z=h.w.E[:, iφ], c=:coolwarm, label=false, markerstrokewidth=0, clims=lims);
for j in axes(w, 2)
    plot!(x, abs2.(w[:, j, iφ]) .+ h.w.E[j, iφ], label=false)
end
plot!(minorgrid=true, xlabel=L"\varphi_x", ylabel="Energy", cbtitle="Energy")