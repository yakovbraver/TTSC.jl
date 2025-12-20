# A script for producing all the figures for https://doi.org/10.1103/PhysRevB.106.144301 (https://arxiv.org/abs/2206.14804)
import TTSC.SineModel as sm
using TTSC.Classical
using Plots, Measures, LaTeXStrings

pythonplot()

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
BLUE2   = colorant"rgb(56, 111, 164)"
GREY2  = colorant"rgb(180, 180, 180)"
BLACK  = colorant"rgb(13, 19, 33)"
CMAP = cgrad(:linear_grey_0_100_c0_n256, rev=true)
hideaxis(x) = "" # use `yformatter=hideaxis` when `yformatter=_->""` randomly doesn't work

l = 1
gₗ = -7640
Vₗ = -2
λₛ = 100; λₗ = 40; ω = 410
s = 2

# swap the Wanniers `w1` and `w2` at phases `ϕ_range`
function swap_wanniers!(w, w1, w2, ϕ_range)
    old1 = w[:, :, w1, ϕ_range]
    w[:, :, w1, ϕ_range] = w[:, :, w2, ϕ_range]
    w[:, :, w2, ϕ_range] = old1
    nothing
end

##########
########## FIG 1
##########

set_defaults(width=2*8.6, height=7.5)

### (a) Floquet spectrum

φₓ = range(0, pi, length=61)
n_cells = 1

h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=34, isperiodic=true)
sm.diagonalise!(h)
H = sm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:time, minband=1)
sm.diagonalise!(H)

figa = plot();
for i in 1:8n_cells
    label, ls, c = (i == 1 ? (L"\beta=1", :solid, BROWN) : i == 3 ? (L"\beta=2", :dash, BROWN) : ("", :solid, GREY))
    plot!(2φₓ ./ π, H.E[i, :]; label, c, ls, lw=0.8)
end
plot!(xlabel=L"\varphi_t/\pi", ylabel="Quasienergy", title="x", legend=(0.01, 0.5))

### (b) Wannier maps

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

targetlevels = [1, 3]
sm.compute_wanniers!(H; targetlevels)

x = range(0, n_cells*pi, length=50n_cells) # x's for wavefunctions
Ωt = range(0, 2π, length=40s) # time moments for wavefunctions
_, w_c = sm.make_wannierfunctions(H, x, Ωt, 1:length(φₓ))
w = abs2.(w_c)

iϕ = 1
fig1 = heatmap(x ./ π, Ωt ./ π, w[:, :, 1, iϕ]', xlabel=L"x/\pi", ylabel=L"\Omega t/\pi", title=L"|w_{1}(x,t)|^2", c=CMAP, cbar=false)
shadecells_2!(fig1)
fig2 = heatmap(x ./ π, Ωt ./ π, w[:, :, 2, iϕ]', xlabel=L"x/\pi", yformatter=_->"", title=L"|w_{2}(x,t)|^2", c=CMAP, cbar=false)
shadecells_2!(fig2)

# with these widths the two plots come out the same in the final plot (probably because the left one has `ylabel` while the right doesn't)
figb = plot(fig1, fig2, layout=grid(1, 2, widths=[0.42, 0.58]), link=:y, xlims=(0, 1), xticks=0:1, ylims=(0, 2), yticks=0:1:2)
figab = plot(figa, figb, layout=grid(2, 1, heights=[0.7, 0.3]))

### (c) Wannier functions

iϕ = 1
i_x = 16
fig1 = plot(Ωt ./ π, w[i_x, :, 1, iϕ], c=BLACK, xformatter=_->"", ylabel=L"|w_\alpha(x_0,t)|^2", label=L"\alpha=1")
plot!(Ωt ./ π, w[i_x, :, 2, iϕ], c=GREY2, label=L"\alpha=2", title=L"\varphi_t=0", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

iϕ = 16
fig2 = plot(Ωt ./ π, w[i_x, :, 1, iϕ], c=BLACK, xlabel=L"\Omega t/\pi", ylabel=L"|w_\alpha(x_0,t)|^2", label=L"\alpha=1")
plot!(Ωt ./ π, w[i_x, :, 2, iϕ], c=GREY2, label=L"\alpha=2", title=L"\varphi_t=\pi/2", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

figc = plot(fig1, fig2, layout=(2,1), link=:x)

### Maps of Wannier functions

fig1 = heatmap(Ωt ./ π, 2φₓ ./ π, w[i_x, :, 1, :]', xformatter=_->"", ylabel=L"\varphi_t/\pi", title="x", cbartitle=L"|w_1(x_0,t)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1], c=GREEN, label=L"\gamma=1", alpha=0.3)
vspan!([1, 2], xformatter=_->"", c=RED, label=L"\gamma=2", alpha=0.3, widen=false, xlims=(0, 2))
fig2 = heatmap(Ωt ./ π, 2φₓ ./ π, w[i_x, :, 2, :]', xlabel=L"\Omega t/\pi", ylabel=L"\varphi_t/\pi", cbartitle=L"|w_2(x_0,t)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1], c=GREEN, label=false, alpha=0.3)
vspan!([1, 2], c=RED, label=false, alpha=0.3, widen=false, xlims=(0, 2))
figd = plot(fig1, fig2, layout=(2,1), link=:x)

plot(figab, figc, figd, layout=(1,3))
savefig("fig1.pdf")

##########
########## FIG 2
##########

set_defaults(width=2*8.6, height=7.5)

### (a) Floquet spectrum

φₓ = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
n_cells = 2

h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=34, isperiodic=true)
sm.diagonalise!(h)
H = sm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:space, minband=1)
sm.diagonalise!(H)

figa = plot();
for i in 1:8n_cells
    label, ls, c = (i == 1 ? (L"j=1,\beta=1", :solid, BROWN) : i == 2 ? (L"j=2,\beta=1", :dash, BROWN) : ("", :solid, GREY))
    plot!(φₓ ./ π, H.E[i, :]; label, c, ls, lw=0.8, xlims=(0, 1))
end
plot!(xlabel=L"\varphi_x/\pi", ylabel="Quasienergy (recoil units)", title="x", legend=(0.01, 0.1), xtick=0:0.25:1, bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))
lens!([0.225, 0.275], [-5696.41, -5696.25], inset = (1, bbox(0.35, 0.25, 0.5, 0.25)), lw=0.5, c=:black)

### (b) Wannier functions

targetlevels = [1, 2]
sm.compute_wanniers!(H; targetlevels)

x = range(0, n_cells*pi, length=1000n_cells) # x's for wavefunctions
Ωt = range(0, 2π, length=40s) # time moments for wavefunctions: 𝜔𝑡/𝑠 ∈ [0; 2π]
_, w_c = sm.make_wannierfunctions(H, x, Ωt, 1:length(φₓ))
w = abs2.(w_c)

i_t = 21

swap_wanniers!(w, 1, 2, 31:lastindex(φₓ))

iϕ = 1
fig1 = plot(x ./ π, w[:, i_t, 1, iϕ], c=BLACK, xformatter=_->"", ylabel=L"|w_{i,1}(x,t_0)|^2", label=L"i=1", lw=0.5)
plot!(x ./ π, w[:, i_t, 2, iϕ], c=GREY2, label=L"i=2", title=L"\varphi_x=0", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3), lw=0.5)

iϕ = 16
fig2 = plot(x ./ π, w[:, i_t, 1, iϕ], c=BLACK, xlabel=L"x/\pi", ylabel=L"|w_{i,1}(x,t_0)|^2", label=L"i=1", lw=0.5, yticks=0:1:2)
plot!(x ./ π, w[:, i_t, 2, iϕ], c=GREY2, label=L"i=2", title=L"\varphi_x=\pi/4", bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3), lw=0.5)

figb = plot(fig1, fig2, layout=(2,1), link=:x)

### (c) Maps of Wannier functions

fig1 = heatmap(x ./ π, φₓ ./ π, w[:, i_t, 1, :]', ylabel=L"\varphi_x/\pi", title="x", cbartitle=L"|w_{1,1}(x,t_0)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1/4-0.02, 5/4+0.02, 7/4-0.02, 7/4+0.02, 2], c=GREEN, label=L"k=1", alpha=0.3)
vspan!([1/4+0.02, 3/4-0.02, 3/4+0.02, 5/4-0.02], xformatter=_->"", c=BLUE, label=L"k=2", alpha=0.3, widen=false, xlims=(0, 2))
fig2 = heatmap(x ./ π, φₓ ./ π, w[:, i_t, 2, :]', xlabel=L"x/\pi", ylabel=L"\varphi_x/\pi", cbartitle=L"|w_{2,1}(x,t_0)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 1/4-0.02, 5/4+0.02, 7/4-0.02, 7/4+0.02, 2], c=GREEN, label=false, alpha=0.3)
vspan!([1/4+0.02, 3/4-0.02, 3/4+0.02, 5/4-0.02], c=BLUE, label=false, alpha=0.3, widen=false, xlims=(0, 2))

figc = plot(fig1, fig2, link=:x, layout=(2,1))
plot(figa, figb, figc, layout=(1,3))

savefig("fig2.pdf")

##########
########## FIG 3
##########

set_defaults(width=2*8.6, height=15)

### (a) Floquet spectrum

φₓ = [range(0, 0.7, length=10); range(0.75, 0.85, length=15); range(0.9, 2.2, length=10); range(2.3, 2.4, length=15); range(2.4, pi, length=10)]
n_cells = 2

h = sm.UnperturbedHamiltonian(n_cells; M=1/2, gₗ, Vₗ, φₓ, maxband=34, isperiodic=true)
sm.diagonalise!(h)
H = sm.FloquetHamiltonian(h; s, λₛ, λₗ, ω, pumptype=:spacetime, minband=1)
sm.diagonalise!(H)

fig1 = plot();
for i in 1:8n_cells
    label, ls, c, = i == 1 ? (L"j=1,\beta=1", :solid, BROWN) : 
        i == 2 ? (L"j=2,\beta=1", :dash, BROWN) :
        i == 5 ? (L"j=1,\beta=2", :dot, BROWN) :
        i == 6 ? (L"j=2,\beta=2", :dashdot, BROWN) : ("", :solid, GREY)
    plot!(2φₓ ./ π, H.E[i, :]; label, c, ls, xlims=(0, 2), lw=0.8)
end
plot!(xlabel=L"\varphi_t=2\varphi_x\ (\pi"*" rad)", ylabel="Quasienergy (recoil units)", title="x", legend=(0.01, 0.35), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))
# lens!([0.475, 0.525], [-5695.9, -5695.8], inset = (1, bbox(0.25, 0.25, 0.5, 0.25)), lw=0.5, c=:black)

### (b)-(f) Maps of Wannier functions

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

function four_wanniers(w, iϕ, ϕ_str, addlabel=false)
    figa = heatmap(x ./ π, Ωt ./ π, w[:, :, 1, iϕ]', xformatter=_->"", ylabel=L"\Omega t/\pi", title=L"|w_{1,1}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figa; addlabel)
    figb = heatmap(x ./ π, Ωt ./ π, w[:, :, 2, iϕ]', xformatter=_->"", yformatter=_->"", title=L"|w_{1,2}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figb; addlabel=false)
    figc = heatmap(x ./ π, Ωt ./ π, w[:, :, 3, iϕ]', xlabel=L"x/\pi", ylabel=L"\Omega t/\pi", title=L"|w_{2,1}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figc; addlabel=false)
    figd = heatmap(x ./ π, Ωt ./ π, w[:, :, 4, iϕ]', xlabel=L"x/\pi", yformatter=_->"", title=L"|w_{2,2}(x,t)|^2", c=CMAP, cbar=false)
    shadecells!(figd; addlabel=false)
    return plot(figa, figb, figc, figd, layout=(2, 2), link=:both, plot_title=L"\varphi_t=2\varphi_x="*ϕ_str, xaxis=((0, 2), 0:1:2), yaxis=((0, 2), 0:1:2))
end

targetlevels = [1, 2, 5, 6]
sm.compute_wanniers!(H; targetlevels)

x = range(0, n_cells*pi, length=50n_cells)
Ωt = range(0, 2π, length=40s)
_, w_c = sm.make_wannierfunctions(H, x, Ωt, 1:length(φₓ))
w = abs2.(w_c)

swap_wanniers!(w, 1, 2, 16)
swap_wanniers!(w, 3, 4, 16)
swap_wanniers!(w, 1, 3, 17:lastindex(φₓ))
swap_wanniers!(w, 2, 4, 17:lastindex(φₓ))

figb = four_wanniers(w, 1, "0")
figc = four_wanniers(w, 16, L"\pi/2")
figd = four_wanniers(w, 31, L"\pi")
fige = four_wanniers(w, 44, L"3\pi/2")
figf = four_wanniers(w, 60, L"2\pi")

# if widths are not specified, then the first column is 1.23 times wider than the others. 
plot(figa, figb, figc, figd, fige, figf, layout=grid(2, 3, widths=[0.285, 0.35, 0.35, 0.285, 0.35, 0.35]), plot_title="")
savefig("fig3.pdf")

##########
########## FIG S1
##########

set_defaults(width=8.6, height=4.3)

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
H_classical = ClassicalHamiltonian(𝐻₀, 𝐻, params, s, min_pos=(1.5, 2), max_pos=(2, 2.5))

Iₛ, M, coeffs = compute_parameters(H_classical, Function[𝑄ₛ, 𝑄ₗ], [2s, s])

Aₛ = abs(coeffs[1]); χₛ = angle(coeffs[1])
Aₗ = abs(coeffs[2]); χₗ = angle(coeffs[2])

### Make a plot of the motion in the (𝐼, ϑ) phase-space in the secular approximation

I_min = 19.5; I_max = 28
ϑ = range(0, 2π, length=100)
I = range(I_min, I_max, length=50)
E = Matrix{Float64}(undef, length(ϑ), length(I))
h₀ = H_classical.𝐸(Iₛ) - ω/s*Iₛ
for i in eachindex(I), t in eachindex(ϑ)
    E[t, i] = h₀ + (I[i]-Iₛ)^2/2M + λₛ*Aₛ*cos(2s*ϑ[t] + χₛ) + λₗ*Aₗ*cos(s*ϑ[t] + χₗ - π/2)
end
figa = contour(ϑ ./ π, I, E', xlabel=L"\Theta/\pi", ylabel=L"I", color=[GREY], minorgrid=true,
               levels=[range(-5725, -5610, length=10); range(-5600, -5575, length=10)], colorbar=false, lw=0.5, xlims=(0, 2))
hline!([Iₛ], c=YELLOW, label=false, lw=0.5)

### Make an "exact" plot of the motion in the (𝐼, ϑ) phase-space

figb = plot();
for (I_min, χ₀) in zip([20, 22.5, 23.5], [0, 0.75, -0.75])
    for i in I_min:0.5:I_max
        I, Θ = compute_IΘ(H_classical, i; n_T=100, χ₀)
        scatter!(Θ ./ π, I, xlabel=L"\theta, rad", markerstrokewidth=0, markeralpha=0.6, label=false, minorgrid=true, c=GREY, markersize=1)
    end
end
plot!(xlabel=L"\Theta/\pi", xlims=(0, 2), yformatter=hideaxis)

plot(figa, figb, link=:y, ylims=(19.5, Inf))
savefig("figS1.pdf")

##########
########## FIG S2
##########

set_defaults(width=2*8.6, height=8.6)

### (a) Quasiclassical energy spectrum

φₜ = range(0, 2π, length=61)
n_cells = s
gₗ = -2λₛ*Aₛ
Vₗ = 2λₗ*Aₗ

h = sm.UnperturbedHamiltonian(n_cells; M, gₗ, Vₗ, φₓ=-φₜ/2, maxband=2, isperiodic=true)
sm.diagonalise!(h)
h.E .+= -(gₗ + Vₗ)/2 + H_classical.𝐸(Iₛ) - ω/s*Iₛ

figa = plot();
for i in 1:2n_cells
    label, ls, c = (i == 1 ? (L"\beta=1", :solid, BROWN) : i == 2 ? (L"\beta=2", :dash, BROWN) : ("", :solid, GREY))
    plot!(φₜ ./ π, h.E[i, :]; label, c, ls, lw=0.8)
end
plot!(xlabel=L"\varphi_t/\pi", ylabel=L"E_{\rm eff}"*" (recoil units)", legend=(0.02, 0.59))

### (b) Quasiclassical Wannier functions

sm.compute_wanniers!(h; targetband=1, mixsubbands=false)
θ = range(0, s*π, length=40s)
_, w = sm.make_wannierfunctions(h, θ, 1:length(φₜ))

𝑈(iϕ) = @. -λₛ*Aₛ*cos(2s*θ) + λₗ*Aₗ*cos(s*θ - φₜ[iϕ]) + H_classical.𝐸(Iₛ) - ω/s*Iₛ

iϕ = 1
fig1 = plot(θ ./ π, 𝑈(iϕ), label=false, c=BLUE2, lw=2)
plot!(θ ./ π, 4abs2.(w[:, 1, iϕ]) .+ h.w.E[1, iϕ], c=BLACK, label=L"|w_1|^2")
plot!(θ ./ π, 4abs2.(w[:, 2, iϕ]) .+ h.w.E[2, iϕ], c=GREY2, label=L"|w_2|^2", ylims=(-5610, -5575))
plot!(xformatter=_->"", ylabel=L"E_{\rm eff}"*" (recoil units)", title=L"\varphi_t=0", legend=(0.01, 0.001), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

iϕ = 16
fig2 = plot(θ ./ π, 𝑈(iϕ), label=false, c=BLUE2, lw=2)
plot!(θ ./ π, 4abs2.(w[:, 1, iϕ]) .+ h.w.E[1, iϕ], c=BLACK, label=L"|w_1|^2")
plot!(θ ./ π, 4abs2.(w[:, 2, iϕ]) .+ h.w.E[2, iϕ], c=GREY2, label=L"|w_2|^2", ylims=(-5610, -5575))
plot!(xlabel=L"\Theta/\pi", ylabel=L"E_{\rm eff}"*" (recoil units)", title=L"\varphi_t=\pi/2", legend=(0.01, 0.001), bgcolorlegend=RGBA(1, 1, 1, 0.3), fgcolorlegend=RGBA(0, 0, 0, 0.3))

figb = plot(fig1, fig2, layout=(2,1), link=:x)

### (c) Maps of quasiclassical Wannier functions

fig1 = heatmap(θ ./ π, φₜ ./ π, abs2.(w[:, 1, :])', xformatter=_->"", ylabel=L"\varphi_t/\pi", title="x", cbartitle=L"|w_1(\Theta)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 0.5, 1.5, 2], c=GREEN, label=L"\gamma=1", alpha=0.3)
vspan!([0.5, 1.5], xformatter=_->"", c=BLUE, label=L"\gamma=2", alpha=0.3, widen=false, xlims=(0, 2), legend=(0.4, 0.2))
fig2 = heatmap(θ ./ π, φₜ ./ π, abs2.(w[:, 2, :])', xlabel=L"\Theta/\pi", ylabel=L"\varphi_t/\pi", cbartitle=L"|w_2(\Theta)|^2", c=CMAP, rightmargin=-10mm)
vspan!([0, 0.5, 1.5, 2], c=GREEN, label=false, alpha=0.3)
vspan!([0.5, 1.5], c=BLUE, label=false, alpha=0.3, widen=false, xlims=(0, 2))
figc = plot(fig1, fig2, layout=(2,1), link=:x)

plot(figa, figb, figc, layout=(1,3))
savefig("figS2.pdf")