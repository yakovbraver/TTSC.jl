using SpecialFunctions: gamma
using Combinatorics: factorial
import Roots
import Optim
import Polynomials
using QuadGK: quadgk

"A type representing spatial Hamiltonian (S2)."
struct SpatialHamiltonian
    l::Int
    g_l::Float64
    V_L::Float64
end

"Construct a spatial Hamiltonian object."
function SpatialHamiltonian(; l::Integer, g::Real, V_L::Real)
    g_l = 2g*factorial(l) / sqrt(π) / gamma(l + 0.5)
    SpatialHamiltonian(l, g_l, V_L)
end

"Periodic potential of Hamiltonian (S2) as a function of position `x`."
function 𝑈(H::SpatialHamiltonian, x::Real)
    H.g_l * cos(2x)^2H.l + H.V_L * cos(x)^2
end

"Momentum of a particle described by Hamiltonian (S2) as a function of position `x` for the given value of energy `E`."
function 𝑝(H::SpatialHamiltonian, E::Real, x::Real)
    p = E - H.g_l * cos(2x)^2H.l - H.V_L * cos(x)^2
    p < 0 ? 0 : sqrt(p) # a safeguard for the case when `x` is slightly outside of the accessible region of oscillations
end

"""
Calculate the action variables for the energies in the first potential minimum to the right of the origin.
Return a tuple of (actions, energies).
"""
function H_of_I(H::SpatialHamiltonian)
    # find position and value of the first potential minimum to the right of the origin
    lo = 0.8; hi = 1.2
    result = Optim.optimize(x -> 𝑈(H, x), lo, hi, Optim.Brent())
    E_min = Optim.minimum(result)
    x_min = Optim.minimizer(result)
    
    # find the value of the first potential maximum to the right of the origin
    lo = 1.3; hi = 1.7
    result = Optim.optimize(x -> -𝑈(H, x), lo, hi, Optim.Brent())
    E_max = -Optim.minimum(result)
    
    n_E = 50 # number of energies (and actions) to save
    I = Vector{Float64}(undef, n_E) # for storing values of the action variable
    E = range(E_min+1e-4, E_max-1e-4, length=n_E) # energies inside the potential "well"

    x_max = x_min # initialise `x_max` -- the second turning point
    for (i, e) in enumerate(E)
        # find the turning points of the oscillations at the given energy `e`
        x_min = Roots.find_zero(x -> 𝑈(H, x) - e, x_min-.05, Roots.Order1(), atol=1e-5)
        x_max = Roots.find_zero(x -> 𝑈(H, x) - e, x_max+.05, Roots.Order1(), atol=1e-5)
        # calculate ∫𝑝d𝑥 for a half-period; the second half is the same hence no division by 2
        I[i] = quadgk(x -> 𝑝(H, e, x), x_min, x_max, rtol=1e-4)[1] / π # `[1]` contains the integral, `[2]` contains error
    end
    I, E
end

"""
Return the first and second derivatives of `y(x)` as a tuple `(y′, y″)`. Use window half-size `ws`.
The values in `y′` and `y″` will correspond to the vector `x`; first and last `ws` values will contain zeros.
"""
function d_and_d²(x, y; ws=2)
    N = length(y)
    y′ = zeros(N)
	y″ = zeros(N)
	for i in 1+ws:N-ws
		f = Polynomials.fit(x[i-ws:i+ws], y[i-ws:i+ws], 2)
		y′[i] = Polynomials.derivative(f, 1)(x[i])
		y″[i] = Polynomials.derivative(f, 2)(x[i])
	end
	y′, y″
end

h0 = SpatialHamiltonian(l=2, g=4.0, V_L=1.0)
I, E = H_of_I(h0)

ws = 2 # window half-size for taking derivatives
E′, E″ = d_and_d²(I, E)

figs = [plot() for _ in 1:3];
figs[1] = plot(I, E, xlabel=L"I", ylabel=L"E");
figs[2] = plot(I[ws+1:end-ws], E′[ws+1:end-ws], xlabel=L"I", ylabel=L"dE/dI", xlims=(I[1], I[end]));
figs[3] = plot(I[ws+1:end-ws], E″[ws+1:end-ws], xlabel=L"I", ylabel=L"dE^2/d^2I", xlims=(I[1], I[end]));
plot(figs..., layout=grid(3,1))