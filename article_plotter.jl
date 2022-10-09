using Plots, Measures, LaTeXStrings

plotlyjs()
pyplot()

"Set plotting defaults and initialise the canvas size with the given `width` and `height` (in cm)"
function set_defaults(;width, height)
    cm2px = 39
    theme(:default, size=(width*cm2px, height*cm2px))
    default(legendfontsize=8, tickfontsize=8, titlefontsize=9, labelfontsize=8, colorbar_titlefontsize=8, plot_titlefontsize=9, xwiden=false, frame=:box)
end

YELLOW = colorant"rgb(219, 173, 106)"
GREY   = colorant"rgb(105, 116, 124)"
RED    = colorant"rgb(237, 71, 74)"
GREEN  = colorant"rgb(132, 221, 99)"
BROWN  = colorant"rgb(78, 1, 16)"   # 224, 186, 215
BLUE   = colorant"rgb(54, 201, 198)"
GREY2  = colorant"rgb(180, 180, 180)"
BLACK  = colorant"rgb(13, 19, 33)"
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)

include("SpacetimeHamiltonian.jl")
include("bandsolvers.jl")

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
λₗAₗ = λₗ*Aₗ; λₛAₛ = λₛ*Aₛ

########## FIG 1

set_defaults(width=2*8.6, height=7.5)

### Floquet spectrum

phases = range(0, pi, length=61)
n_cells = 1
n_max = 34
n_target = 1
x = range(0, n_cells*pi, length=50n_cells) # x's for wavefunctions
ωts = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
e, E, pos_lo, pos_hi, ε_lo, ε_hi, wf_lo, wf_hi, u_lo, u_hi = compute_floquet_wannier_centres(;N=n_cells, n_target, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, coords=x, ωts, pumptype=:time)

fig1 = plot();
for i in 1:8n_cells
    label, ls, c = (i == 1 ? (L"\beta=1", :solid, BROWN) : i == 3 ? (L"\beta=2", :dash, BROWN) : ("", :solid, GREY))
    plot!(2phases ./ π, E[i, :]; label, c, ls, lw=0.8)
end
plot!(xlabel=L"\varphi_t/\pi", ylabel="Quasienergy", title="x", legend=(0.01, 0.5))

### Wannier maps

function shadecells_2!(fig)
    fillalpha = 0.3; alpha = 0; label=false
    plot!(fig, [0, 1/4-0.02], [2, 2]; fillrange=1, fillalpha, alpha, c=GREEN, label)
    plot!(fig, [0, 1/4-0.02], [1, 1]; fillrange=0, fillalpha, alpha, c=RED, label)
    
    plot!(fig, [1/4+0.02, 2/4], [2, 2]; fillrange=1, fillalpha, alpha, c=RED, label)
    plot!(fig, [1/4+0.02, 2/4], [1, 1]; fillrange=0, fillalpha, alpha, c=GREEN, label)
    plot!(fig, [2/4, 3/4-0.02], [2, 2]; fillrange=1, fillalpha, alpha, c=GREEN, label)
    plot!(fig, [2/4, 3/4-0.02], [1, 1]; fillrange=0, fillalpha, alpha, c=RED, label)
    
    plot!(fig, [3/4+0.02, 4/4], [2, 2]; fillrange=1, fillalpha, alpha, c=RED, label)
    plot!(fig, [3/4+0.02, 4/4], [1, 1]; fillrange=0, fillalpha, alpha, c=GREEN, label)
end

i_ϕ = 1
figa = heatmap(x ./ π, ωts ./ π, wf_hi[:, 2, :, i_ϕ]', xlabel=L"x/\pi", ylabel=L"\Omega t/\pi", title=L"|w_{1}(x,t)|^2", c=CMAP, cbar=false)
shadecells_2!(figa)
figb = heatmap(x ./ π, ωts ./ π, wf_hi[:, 1, :, i_ϕ]', xlabel=L"x/\pi", yformatter=_->"", title=L"|w_{2}(x,t)|^2", c=CMAP, cbar=false)
shadecells_2!(figb)

# with these widths the two plots come out the same (probably because the left on has ylabel while the right doesn't)
fig2n = plot(figa, figb, layout=grid(1, 2, widths=[0.42, 0.58]), link=:y, xlims=(0, 1), xticks=0:1, ylims=(0, 2), yticks=0:1:2)

# floquet modes
figa = heatmap(x ./ π, ωts ./ π, u_hi[:, 2, :, i_ϕ]', xlabel=L"x/\pi", ylabel=L"\omega t/s\ (\pi"*" rad)", title=L"|u_{1}(x,t)|^2", c=CMAP, cbar=false)
figb = heatmap(x ./ π, ωts ./ π, u_hi[:, 1, :, i_ϕ]', xlabel=L"x/\pi", ylabel=L"\omega t/s\ (\pi"*" rad)", title=L"|u_{2}(x,t)|^2", c=CMAP, cbar=false)
figu = plot(figa, figb, layout=grid(1, 2, widths=[0.5, 0.5]), link=:y, xlims=(0, 1), xticks=0:1, ylims=(0, 2), yticks=0:1:2)
savefig("floqs.pdf")
### Wannier functions

i_ϕ = 1
i_x = 16
fig2 = plot(ωts ./ π, wf_hi[i_x, 2, :, i_ϕ], c=BLACK, xformatter=_->"", ylabel=L"|w_\alpha(x_0,t)|^2", label=L"\alpha=1")
plot!(ωts ./ π, wf_hi[i_x, 1, :, i_ϕ], c=GREY2, label=L"\alpha=2", title=L"\varphi_t=0", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

i_ϕ = 16
i_x = 16
fig3 = plot(ωts ./ π, wf_hi[i_x, 2, :, i_ϕ], c=BLACK, xlabel=L"\Omega t/\pi", ylabel=L"|w_\alpha(x_0,t)|^2", label=L"\alpha=1")
plot!(ωts ./ π, wf_hi[i_x, 1, :, i_ϕ], c=GREY2, label=L"\alpha=2", title=L"\varphi_t=\pi/2", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

fig23 = plot(fig2, fig3, layout=(2,1), link=:x)


### Maps of Wannier functions

fig4 = heatmap(ωts ./ π, 2phases ./ π, wf_hi[i_x, 2, :, :]', xformatter=_->"", ylabel=L"\varphi_t/\pi", title="x", cbartitle=L"|w_1(x_0,t)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1], c=GREEN, label=L"\gamma=1", alpha=0.3)
vspan!([1, 2], xformatter=_->"", c=RED, label=L"\gamma=2", alpha=0.3, widen=false, xlims=(0, 2))
fig5 = heatmap(ωts ./ π, 2phases ./ π, wf_hi[i_x, 1, :, :]', xlabel=L"\Omega t/\pi", ylabel=L"\varphi_t/\pi", cbartitle=L"|w_2(x_0,t)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1], c=GREEN, label=false, alpha=0.3)
vspan!([1, 2], c=RED, label=false, alpha=0.3, widen=false, xlims=(0, 2))
fig45 = plot(fig4, fig5, layout=(2,1), link=:x)



fig12 = plot(fig1, fig2n, layout=grid(2, 1, heights=[0.7, 0.3]))

plot(fig12, fig23, fig45, layout=(1,3))
savefig("fig1.pdf")

########## FIG 2

set_defaults(width=2*8.6, height=7.5)

### Floquet spectrum

phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
n_cells = 2
n_max = 34
n_target = 1
x = range(0, n_cells*pi, length=1000n_cells) # x's for wavefunctions
ωts = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
e, E, pos_lo, pos_hi, ε_lo, ε_hi, wf_lo, wf_hi = compute_floquet_wannier_centres(;N=n_cells, n_target, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, coords=x, ωts, mix_time_cells=false, pumptype=:space)

fig1 = plot();
for i in 1:8n_cells
    label, ls, c = (i == 1 ? (L"j=1,\beta=1", :solid, BROWN) : i == 2 ? (L"j=2,\beta=1", :dash, BROWN) : ("", :solid, GREY))
    plot!(phases ./ π, E[i, :]; label, c, ls, lw=0.8, xlims=(0, 1))
end
plot!(xlabel=L"\varphi_x/\pi", ylabel="Quasienergy (recoil units)", title="x", legend=(0.01, 0.1), xtick=0:0.25:1, bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))
lens!([0.225, 0.275], [-5696.41, -5696.25], inset = (1, bbox(0.35, 0.25, 0.5, 0.25)), lw=0.5, c=:black)

### Wannier functions

i_t = 21

# swap the Wanniers 1 and 2 at phases 31:end
old1 = wf_hi[:, 1, i_t, 31:end]
wf_hi[:, 1, i_t, 31:end] = wf_hi[:, 2, i_t, 31:end]
wf_hi[:, 2, i_t, 31:end] = old1

i_ϕ = 1
fig2 = plot(x ./ π, wf_hi[:, 1, i_t, i_ϕ], c=YELLOW, xformatter=_->"", ylabel=L"|w_{i,1}(x,t_0)|^2", label=L"i=1", lw=0.5)
plot!(x ./ π, wf_hi[:, 2, i_t, i_ϕ], c=RED, label=L"i=2", title=L"\varphi_x=0", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3), lw=0.5)

i_ϕ = 16
fig3 = plot(x ./ π, wf_hi[:, 1, i_t, i_ϕ], c=YELLOW, xlabel=L"x/\pi", ylabel=L"|w_{i,1}(x,t_0)|^2", label=L"i=1", lw=0.5, yticks=0:1:2)
plot!(x ./ π, wf_hi[:, 2, i_t, i_ϕ], c=RED, label=L"i=2", title=L"\varphi_x=\pi/4", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3), lw=0.5)

fig23 = plot(fig2, fig3, layout=(2,1), link=:x)

### Maps of Wannier functions

fig4 = heatmap(x ./ π, phases ./ π, wf_hi[:, 1, i_t, :]', ylabel=L"\varphi_x/\pi", title="x", cbartitle=L"|w_{1,1}(x,t_0)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1/4-0.02, 5/4+0.02, 7/4-0.02, 7/4+0.02, 2], c=GREEN, label=L"k=1", alpha=0.3)
vspan!([1/4+0.02, 3/4-0.02, 3/4+0.02, 5/4-0.02], xformatter=_->"", c=BLUE, label=L"k=2", alpha=0.3, widen=false, xlims=(0, 2))
fig5 = heatmap(x ./ π, phases ./ π, wf_hi[:, 2, i_t, :]', xlabel=L"x/\pi", ylabel=L"\varphi_x/\pi", cbartitle=L"|w_{2,1}(x,t_0)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1/4-0.02, 5/4+0.02, 7/4-0.02, 7/4+0.02, 2], c=GREEN, label=false, alpha=0.3)
vspan!([1/4+0.02, 3/4-0.02, 3/4+0.02, 5/4-0.02], c=BLUE, label=false, alpha=0.3, widen=false, xlims=(0, 2))

fig45 = plot(fig4, fig5, link=:x, layout=(2,1))

plot(fig1, fig23, fig45, layout=(1,3))

savefig("fig2.pdf")

########## FIG 3

set_defaults(width=2*8.6, height=7.5)

### Floquet spectrum

phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
n_cells = 2
n_max = 34
n_target = 1
x = range(0, n_cells*pi, length=50n_cells) # x's for wavefunctions
ωts = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
e, E, pos_lo, pos_hi, ε_lo, ε_hi, wf_lo, wf_hi, u_lo, u_hi = compute_floquet_wannier_centres(;N=n_cells, n_target, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, coords=x, ωts, pumptype=:spacetime)

fig1 = plot();
for i in 1:8n_cells
    label, ls, c, = i == 1 ? (L"j=1,\beta=1", :solid, BROWN) : 
        i == 2 ? (L"j=2,\beta=1", :dash, BROWN) :
        i == 5 ? (L"j=1,\beta=2", :dot, BROWN) :
        i == 6 ? (L"j=2,\beta=2", :dashdot, BROWN) : ("", :solid, GREY)
    plot!(2phases ./ π, E[i, :]; label, c, ls, xlims=(0, 2), lw=0.8)
end
plot!(xlabel=L"\varphi_t=2\varphi_x\ (\pi"*" rad)", ylabel="Quasienergy (recoil units)", title="x", legend=(0.01, 0.35), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))
# lens!([0.475, 0.525], [-5695.9, -5695.8], inset = (1, bbox(0.25, 0.25, 0.5, 0.25)), lw=0.5, c=:black)

### Maps of Wannier functions

# swap the Wanniers `w1` and `w2` at phases `ϕ_range`
function swap_wanniers(w1, w2, ϕ_range)
    old1 = wf_hi[:, w1, :, ϕ_range]
    wf_hi[:, w1, :, ϕ_range] = wf_hi[:, w2, :, ϕ_range]
    wf_hi[:, w2, :, ϕ_range] = old1
    nothing
end

swap_wanniers(1, 2, 16)
swap_wanniers(3, 4, 16)
swap_wanniers(1, 3, 17:lastindex(phases))
swap_wanniers(2, 4, 17:lastindex(phases))

function shadecells!(fig; addlabel)
    fillalpha = 0.3; alpha = 0; 
    label = addlabel ? L"k=1,\gamma=1" : false
    plot!(fig, [0, 1/4-0.02], [2, 2]; fillrange=1, fillalpha, alpha, c=GREEN, label)
    plot!(fig, [7/4+0.02, 2], [1, 1]; fillrange=0, fillalpha, alpha, c=GREEN, label=false)
    label = addlabel ? L"k=1,\gamma=2" : false
    plot!(fig, [0, 1/4-0.02], [1, 1]; fillrange=0, fillalpha, alpha, c=RED, label)
    plot!(fig, [7/4+0.02, 2], [2, 2]; fillrange=1, fillalpha, alpha, c=RED, label=false)
    
    plot!(fig, [5/4+0.02, 6/4], [1, 1]; fillrange=0, fillalpha, alpha, c=GREEN, label=false)
    plot!(fig, [6/4, 7/4-0.02], [2, 2]; fillrange=1, fillalpha, alpha, c=GREEN, label=false)
    plot!(fig, [5/4+0.02, 6/4], [2, 2]; fillrange=1, fillalpha, alpha, c=RED, label=false)
    plot!(fig, [6/4, 7/4-0.02], [1, 1]; fillrange=0, fillalpha, alpha, c=RED, label=false)
    
    label = addlabel ? L"k=2,\gamma=1" : false
    plot!(fig, [1/4+0.02, 2/4], [1, 1]; fillrange=0, fillalpha, alpha, c=YELLOW, label)
    plot!(fig, [2/4, 3/4-0.02], [2, 2]; fillrange=1, fillalpha, alpha, c=YELLOW, label=false)
    label = addlabel ? L"k=2,\gamma=2" : false
    plot!(fig, [1/4+0.02, 2/4], [2, 2]; fillrange=1, fillalpha, alpha, c=BLUE, label)
    plot!(fig, [2/4, 3/4-0.02], [1, 1]; fillrange=0, fillalpha, alpha, c=BLUE, label=false)
    
    plot!(fig, [3/4+0.02, 4/4], [2, 2]; fillrange=1, fillalpha, alpha, c=BLUE, label=false)
    plot!(fig, [4/4, 5/4-0.02], [1, 1]; fillrange=0, fillalpha, alpha, c=BLUE, label=false)
    plot!(fig, [3/4+0.02, 4/4], [1, 1]; fillrange=0, fillalpha, alpha, c=YELLOW, label=false)
    plot!(fig, [4/4, 5/4-0.02], [2, 2]; fillrange=1, fillalpha, alpha, c=YELLOW, label=false)
end

function four_wanniers(i_ϕ, ϕ_str, addlabel=false)
    figa = heatmap(x ./ π, ωts ./ π, wf_hi[:, 1, :, i_ϕ]', xformatter=_->"", ylabel=L"\Omega t/\pi", title=L"|w_{1,1}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figa; addlabel)
    figb = heatmap(x ./ π, ωts ./ π, wf_hi[:, 2, :, i_ϕ]', xformatter=_->"", yformatter=_->"", title=L"|w_{1,2}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figb; addlabel=false)
    figc = heatmap(x ./ π, ωts ./ π, wf_hi[:, 3, :, i_ϕ]', xlabel=L"x/\pi", ylabel=L"\Omega t/\pi", title=L"|w_{2,1}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figc; addlabel=false)
    figd = heatmap(x ./ π, ωts ./ π, wf_hi[:, 4, :, i_ϕ]', xlabel=L"x/\pi", yformatter=_->"", title=L"|w_{2,2}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figd; addlabel=false)
    return plot(figa, figb, figc, figd, layout=(2, 2), link=:both, plot_title=L"\varphi_t=2\varphi_x="*ϕ_str, xaxis=((0, 2), 0:1:2), yaxis=((0, 2), 0:1:2))
end

fig2 = four_wanniers(1, "0", true)
fig3 = four_wanniers(60, L"2\pi")

# if widths are not specified, then the first column is 1.23 times wider than the others. 
plot(fig1, fig2, fig3, layout=grid(1, 3, widths=[0.285, 0.35, 0.35]), plot_title="")
savefig("fig4.pdf")

########## FIG S1

set_defaults(width=8.6, height=4.3)

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation

function plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, φₜ, Iₛ, s, I_min, I_max)
    ϑ = range(0, 2π, length=100)
    I = range(I_min, I_max, length=50)
    E = Matrix{Float64}(undef, length(ϑ), length(I))
    h₀ = H.𝐸(Iₛ) - ω/s*Iₛ
    for i in eachindex(I), t in eachindex(ϑ)
        E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] + χₛ) + λₗ*Aₗ*cos(s*ϑ[t] + χₗ - φₜ)
    end
    contour(ϑ ./ π, I, E', xlabel=L"\Theta/\pi", ylabel=L"I", color=[GREY], minorgrid=true,
           levels=[range(-5725, -5610, length=10); range(-5600, -5575, length=10)], colorbar=false, lw=0.5, xlims=(0, 2))
    hline!([Iₛ], c=YELLOW, label=false, lw=0.5)
end

I_min = 19.5; I_max = 28
fig1 = plot_isoenergies(; ω, M, λₛ, Aₛ, χₛ, λₗ, Aₗ, χₗ, φₜ=π/2, Iₛ, s, I_min, I_max)

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

fig2 = plot();
for (I_min, χ₀) in zip([20, 22.5, 23.5], [0, 0.75, -0.75])
    for i in I_min:0.5:I_max
        I, Θ = compute_IΘ(H, i; n_T=100, χ₀) # use χ₀ = 0 and ±0.75
        scatter!(Θ ./ π, I, xlabel=L"\theta, rad", markerstrokewidth=0, markeralpha=0.6, label=false, minorgrid=true, c=GREY, markersize=1)
    end
end
plot!(xlabel=L"\Theta/\pi", xlims=(0, 2), yformatter=_->"")

plot(fig1, fig2, link=:y, ylims=(19.5, Inf))
savefig("figS1.pdf")

########## FIG S2

set_defaults(width=2*8.6, height=8.6)

### Quasiclassical energy spectrum

φₜ = range(0, 2π, length=61)
n_cells = s
gₗ = -2λₛ*Aₛ
Vₗ = 2λₗ*Aₗ

h = Bandsolvers.UnperturbedHamiltonian(n_cells; M, gₗ, Vₗ, phases=-φₜ/2, maxband=1, isperiodic=true)
Bandsolvers.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H.𝐸(Iₛ) - ω/s*Iₛ

fig1 = plot();
for i in 1:2n_cells
    label, ls, c = (i == 1 ? (L"\beta=1", :solid, BROWN) : i == 2 ? (L"\beta=2", :dash, BROWN) : ("", :solid, GREY))
    plot!(φₜ ./ π, h.E[i, :]; label, c, ls, lw=0.8)
end
plot!(xlabel=L"\varphi_t/\pi", ylabel=L"E_{\rm eff}"*" (recoil units)", legend=(0.02, 0.59))


### Quasiclassical Wannier functions

Bandsolvers.compute_wanniers!(h, targetband=1)
θ = range(0, s*π, length=40s)
w_lo, _ = Bandsolvers.make_wannierfunctions(h, θ, 1:length(φₜ))

𝑈(i_ϕ) = @. -λₛ*Aₛ*cos(2s*θ) + λₗ*Aₗ*cos(s*θ - φₜ[i_ϕ]) + H.𝐸(Iₛ) - ω/s*Iₛ

i_ϕ = 1
fig2 = plot(θ ./ π, 𝑈(i_ϕ), label=false, c=GREY, lw=2)
plot!(θ ./ π, 4abs2.(w_lo[i_ϕ][1]) .+ h.w.E_lo[i_ϕ][1], c=YELLOW, label=L"|w_1|^2")
plot!(θ ./ π, 4abs2.(w_lo[i_ϕ][2]) .+ h.w.E_lo[i_ϕ][2], c=RED, label=L"|w_2|^2", ylims=(-5610, -5575))
plot!(xformatter=_->"", ylabel=L"E_{\rm eff}"*" (recoil units)", title=L"\varphi_t=0", legend=(0.01, 0.001), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

i_ϕ = 16
fig3 = plot(θ ./ π, 𝑈(i_ϕ), label=false, c=GREY, lw=2)
plot!(θ ./ π, 4abs2.(w_lo[i_ϕ][1]) .+ h.w.E_lo[i_ϕ][1], c=YELLOW, label=L"|w_1|^2")
plot!(θ ./ π, 4abs2.(w_lo[i_ϕ][2]) .+ h.w.E_lo[i_ϕ][2], c=RED, label=L"|w_2|^2", ylims=(-5610, -5575))
plot!(xlabel=L"\Theta/\pi", ylabel=L"E_{\rm eff}"*" (recoil units)", title=L"\varphi_t=\pi/2", legend=(0.01, 0.001), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

fig23 = plot(fig2, fig3, layout=(2,1), link=:x)

### Maps of quasiclassical Wannier functions

w_map = Bandsolvers.make_wanniermap(w_lo, 1) .|> abs2
fig4 = heatmap(θ ./ π, φₜ ./ π, w_map', xformatter=_->"", ylabel=L"\varphi_t/\pi", title="x", cbartitle=L"|w_1(\Theta)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 0.5, 1.5, 2], c=GREEN, label=L"\gamma=1", alpha=0.3)
vspan!([0.5, 1.5], xformatter=_->"", c=BLUE, label=L"\gamma=2", alpha=0.3, widen=false, xlims=(0, 2), legend=(0.4, 0.2))
w_map = Bandsolvers.make_wanniermap(w_lo, 2) .|> abs2
fig5 = heatmap(θ ./ π, φₜ ./ π, w_map', xlabel=L"\Theta/\pi", ylabel=L"\varphi_t/\pi", cbartitle=L"|w_2(\Theta)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 0.5, 1.5, 2], c=GREEN, label=false, alpha=0.3)
vspan!([0.5, 1.5], c=BLUE, label=false, alpha=0.3, widen=false, xlims=(0, 2))
fig45 = plot(fig4, fig5, layout=(2,1), link=:x)

plot(fig1, fig23, fig45, layout=(1,3))
savefig("figS2.pdf")

########## FIG S3

set_defaults(width=2*8.6, height=15)

phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
n_cells = 2
n_max = 34
n_target = 1
x = range(0, n_cells*pi, length=50n_cells) # x's for wavefunctions
ωts = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
e, E, pos_lo, pos_hi, ε_lo, ε_hi, wf_lo, wf_hi = compute_floquet_wannier_centres(;N=n_cells, n_target, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, coords=x, ωts, pumptype=:spacetime)

fig6 = plot();
for i in 1:8n_cells
    label, ls, c, = i == 1 ? (L"j=1,\beta=1", :solid, BROWN) : 
        i == 2 ? (L"j=2,\beta=1", :dash, BROWN) :
        i == 5 ? (L"j=1,\beta=2", :dot, BROWN) :
        i == 6 ? (L"j=2,\beta=2", :dashdot, BROWN) : ("", :solid, GREY)
    plot!(2phases ./ π, E[i, :]; label, c, ls, xlims=(0, 2), lw=0.8)
end
plot!(xlabel=L"\varphi_t=2\varphi_x\ (\pi"*" rad)", ylabel="Quasienergy (recoil units)", title="x", legend=(0.01, 0.35), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))


### Maps of Wannier functions

phases = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
n_cells = 2
n_max = 34
n_target = 1
x = range(0, n_cells*pi, length=50n_cells) # x's for wavefunctions
ωts = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
e, E, pos_lo, pos_hi, ε_lo, ε_hi, wf_lo, wf_hi = compute_floquet_wannier_centres(;N=n_cells, n_target, n_max, phases, s, gₗ, Vₗ, λₗ, λₛ, ω, coords=x, ωts, pumptype=:spacetime)

swap_wanniers(1, 2, 16)
swap_wanniers(3, 4, 16)
swap_wanniers(1, 3, 17:lastindex(phases))
swap_wanniers(2, 4, 17:lastindex(phases))

fig1 = four_wanniers(1, "0")
fig2 = four_wanniers(16, L"\pi/2")
fig3 = four_wanniers(31, L"\pi")
fig4 = four_wanniers(44, L"3\pi/2")
fig5 = four_wanniers(60, L"2\pi")
fig6 = plot()

# if widths are not specified, then the first column is 1.23 times wider than the others. 
# plot(fig1, fig2, fig3, fig4, fig5, fig6, layout=grid(2, 3), plot_title="")
savefig("figS3.pdf")
plot(fig6, fig1, fig2, fig3, fig4, fig5, layout=grid(2, 3, widths=[0.285, 0.35, 0.35, 0.285, 0.35, 0.35]), plot_title="")
savefig("fig4.pdf")