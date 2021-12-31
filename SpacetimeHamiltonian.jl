import Roots
import Optim
import Dierckx
import Elliptic
using QuadGK: quadgk
import DifferentialEquations as DiffEq
using DiffEqPhysics: HamiltonianProblem, DynamicalODEProblem

"A type representing spatial Hamiltonian."
mutable struct SpacetimeHamiltonian
    𝐻₀::Function     # free (unperturbed) Hamiltonian
    𝑈::Function # spatial potential
    left_tp::Tuple
    right_tp::Tuple
    𝐸::Dierckx.Spline1D # energy at the given angle, function 𝐸(𝐼)
    𝐸′::Function    # oscillation frequency at the given angle, function 𝐸′(𝐼)
    𝐸″::Function    # effective mass at the given angle, function 𝐸″(𝐼)
    d𝑥╱d𝑡!::Function
    d𝑝╱d𝑡!::Function
    params::Vector{Float64} # a vector of parameters, will be shared among 𝐻₀, d𝑥╱d𝑡!, and d𝑝╱d𝑡!; the last element should contain external frequency
    s::Int  # resonance number
end

function SpacetimeHamiltonian(𝐻₀::Function, left_tp::Tuple{<:Real, <:Real}, right_tp::Tuple{<:Real, <:Real}, min_pos::Tuple{<:Real, <:Real}, max_pos::Tuple{<:Real, <:Real},
                              d𝑥╱d𝑡!::Function, d𝑝╱d𝑡!::Function, params::AbstractVector, s::Integer)
    𝑈 = x -> 𝐻₀(0.0, x, params)
    𝐸, 𝐸′, 𝐸″ = make_action_functions(𝑈, min_pos, max_pos)
    SpacetimeHamiltonian(𝐻₀, 𝑈, left_tp, right_tp, 𝐸, 𝐸′, 𝐸″, d𝑥╱d𝑡!, d𝑝╱d𝑡!, params, s)
end

"Momentum 𝑝(𝑥) = √[𝐸 - 𝑈(𝑥)] of a particle of energy `E`."
function 𝑝(𝑈::Function, E::Real, x::Real)
    p = E - 𝑈(x)
    p < 0 ? zero(p) : sqrt(p) # a safeguard for the case when `x` is slightly outside of the accessible region of oscillations
end

"""
Construct the functions 𝐸(𝐼), 𝐸′(𝐼), and 𝐸″(𝐼) witing them into `H`.
"""
function make_action_functions(𝑈::Function, min_pos::Tuple{<:Real, <:Real}, max_pos::Tuple{<:Real, <:Real})
    # find position and value of the potential minimum
    result = Optim.optimize(x -> 𝑈(first(x)), min_pos[1], min_pos[2], Optim.Brent())
    E_min = Optim.minimum(result)
    x_min = Optim.minimizer(result)
    
    # find the value of the potential maximum
    result = Optim.optimize(x -> -𝑈(first(x)), max_pos[1], max_pos[2], Optim.Brent())
    E_max = -Optim.minimum(result)
    
    n_E = 100 # number of energies (and actions) to save
    I = Vector{Float64}(undef, n_E) # for storing values of the action variable
    E = range(E_min+1e-4, E_max-1e-4, length=n_E) # energies inside the potential "well"

    x_max = x_min # initialise `x_max` -- the second turning point
    for i in eachindex(E)
        # we find the turning points manually beacuse we need them back for the next iteration
        x_min, x_max = turning_points(𝑈, E[i], x_min - 0.05, x_max + 0.05)
        I[i] = 𝐼(𝑈, E[i], (x_min, x_max))
    end
    
    𝐸  = Dierckx.Spline1D(I, E; k=2)      
    𝐸′ = x -> Dierckx.derivative(𝐸, x; nu=1)
    𝐸″ = x -> Dierckx.derivative(𝐸, x; nu=2)
    return 𝐸, 𝐸′, 𝐸″
end

"""
Return the turning points of motion with the given energy `E`. The initial guesses `a` and `b` may be given as single numbers
or as tuples representing the bracketing interval.
"""
function turning_points(𝑈::Function, E::Real, a::Union{<:Real, Tuple{<:Real, <:Real}}, b::Union{<:Real, Tuple{<:Real, <:Real}})
    x_min = Roots.find_zero(x -> 𝑈(x) - E, a, atol=1e-5)
    x_max = Roots.find_zero(x -> 𝑈(x) - E, b, atol=1e-5)
    return x_min, x_max
end

"""
Return action for the given energy `E` as the integral of momentum over a period of motion.
The turning points should be specified as a tuple `turnpoints`.
"""
function 𝐼(𝑈::Function, E::Real, turnpoints::Tuple{<:Real, <:Real})
    x_min, x_max = turnpoints
    # calculate ∫𝑝d𝑥 for a half-period; the second half is the same hence no division by 2
    return quadgk(x -> 𝑝(𝑈, E, x), x_min, x_max, rtol=1e-4)[1] / π # `[1]` contains the integral, `[2]` contains error
end

"""
Return action for the given energy `E` as the integral of momentum over a period of motion.
The turning points will be determined using the bracketing intervals `H.left_tp` and `H.right_tp`.
"""
function 𝐼(H::SpacetimeHamiltonian, E::Real)
    x_min, x_max = turning_points(H.𝑈, E, H.left_tp, H.right_tp)
    # calculate ∫𝑝d𝑥 for a half-period; the second half is the same hence no division by 2
    return quadgk(x -> 𝑝(H.𝑈, E, x), x_min, x_max, rtol=1e-4)[1] / π # `[1]` contains the integral, `[2]` contains error
end

"Return the action and mass at the working point, as well as the `H.s`th Fourier coefficients for every function in `perturbations`."
function compute_parameters(H::SpacetimeHamiltonian, perturbations::Vector{Function})
    ω = H.params[end]
    Ω = ω / H.s # our choice of the oscillation frequency (of the unperturbed system)
    Iₛ::Float64 = Roots.find_zero(x -> H.𝐸′(x) - Ω, (0, Dierckx.get_knots(H.𝐸)[end]), atol=1e-5) # find which 𝐼ₛ gives the frequency Ω
    E₀::Float64 = H.𝐸(Iₛ)     # energy of the system oscillating at the frequency Ω
    M::Float64 = 1 / H.𝐸″(Iₛ) # "mass" of the system oscillating at the frequency Ω

    # Evolve the unperturbed system for one period 
    T = 2π / Ω
    tspan = (0.0, T) # we use the theoretical value of the period
    # Initial conditions; they may be chosen arbitrary as long as the total energy equals `E₀`
    q₀ = 0.0; p₀ = sqrt(E₀); # We assume that 𝐸 = 𝑝² + 𝑈(𝑥) with 𝑈(0) = 0
    H₀_problem = HamiltonianProblem(H.𝐻₀, p₀, q₀, tspan, H.params)
    dt=2e-4
    # none of RKN solvers worked (https://docs.juliahub.com/DifferentialEquations/UQdwS/6.15.0/solvers/dynamical_solve/)
    sol = DiffEq.solve(H₀_problem, DiffEq.McAte3(); dt) # McAte3 is more accurate than the automatically chosen Tsit5() 

    # calculate 𝑠th Fourier coefficient for every function in `perturbations`
    coeffs = Vector{Float64}(undef, length(perturbations))
    V = Vector{Float64}(undef, length(sol.t)) # for storing perturbation evaluated in the solution points
    for (i, 𝑉) in enumerate(perturbations)
        V .= 𝑉.(sol[2, :], sol[1, :])
        coeffs[i] = fourier_coeff(V, s, dt, T) |> abs
    end
    return Iₛ, M, coeffs
end

"Calculate the `n`th Fourier coefficient of `f`. Simple trapezoid rule is used."
function fourier_coeff(f::AbstractVector, n::Int, dt::AbstractFloat, T::AbstractFloat)
    (sum(f[i] * cispi(2n*(i-1)*dt/T) for i = 2:length(f)-1) + (f[1] + f[end])/2) * dt/T
end

function compute_IΘ(H::SpacetimeHamiltonian, I_target::Real)    
    ω = H.params[end]
    T_external = 2π / ω # period of the external driving
    n_T = 50 # number of periods of the external driving to calculate evolution for
    tspan = (0.0, n_T * T_external)
    
    x₀ = 0.0; p₀ = sqrt(H.𝐸(I_target));
    H_problem = DynamicalODEProblem(H.d𝑥╱d𝑡!, H.d𝑝╱d𝑡!, [x₀], [p₀], tspan, params)
    sol = DiffEq.solve(H_problem, DiffEq.KahanLi8(); dt=2e-4, saveat=T_external)

    x = sol[1, :]
    p = sol[2, :]
    E = map((p, x) -> H.𝐻₀(p, x, H.params), p, x)
    I = map(x -> 𝐼(H, x), E)

    # Find phases from the coordinates
    Θ = similar(I)
    for i in eachindex(Θ)
        T_free = 2π / H.𝐸′(I[i]) # period of the unperturbed motion at action `I[i]`

        tspan = (0.0, T_free)
        # Initial conditions; they may be chosen arbitrary as long as the total energy equals `E₀`
        q₀ = 0.0; p₀ = sqrt(E[i]);
        H₀_problem = HamiltonianProblem(H.𝐻₀, p₀, q₀, tspan, H.params)
        sol = DiffEq.solve(H₀_problem, DiffEq.McAte3(); dt=2e-4)
        # plot(sol) |> display
        
        bracket = x[i] > 0 ? (T_free/2, T_free) : (0.0, T_free/2)
        t = Roots.find_zero(t -> sol(t)[1] - p[i], bracket, Roots.A42(), xrtol=1e-3)
        
        Θ[i] = rem2pi(t / T_free * 2π - pi/2, RoundDown) # `-2π*(i-1)/H.s` is the -ω𝑡/𝑠 term that transforms to the moving frame. We have 𝑡ₙ = 𝑛𝑇, and ω𝑡ₙ = 2π𝑛
    end
    return I, Θ
end