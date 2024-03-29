module Classical

import Roots
import Optim
import Dierckx
using QuadGK: quadgk
import OrdinaryDiffEq as DiffEq
using DiffEqPhysics: HamiltonianProblem

"A type representing a classical Hamiltonian."
mutable struct ClassicalHamiltonian
    𝐻₀::Function    # free (unperturbed) Hamiltonian
    𝐻::Function     # full Hamiltonian, including time-dependent perturbation
    𝑈::Function     # spatial potential
    left_tp::Tuple{Float64, Float64}  # bracketing interval for the left turning point of the free motion
    right_tp::Tuple{Float64, Float64} # bracketing interval for the right turning point of the free motion
    𝐸::Dierckx.Spline1D # energy at the given action, function 𝐸(𝐼)
    𝐸′::Function    # oscillation frequency at the given action, function 𝐸′(𝐼)
    𝐸″::Function    # effective mass at the given action, function 𝐸″(𝐼)
    params::Vector{Float64} # a vector of parameters, will be shared among 𝐻₀ and 𝐻; the last element should contain external frequency
    s::Int  # resonance number
end

"""
Construct a `ClassicalHamiltonian` object. Either known intervals for the turning points have to be provided (`left_tp` and `right_tp`),
or the bracketing intervals for the minimum and the maximum (`min_pos` and `max_pos`) of the spatial potential so that the turning point
intervals are determined automatically. In the latter case, `turnpoint` is required if the potential is not symmetric, see [`turning_point_intervals`](@ref).
"""
function ClassicalHamiltonian(𝐻₀::Function, 𝐻::Function, params::AbstractVector, s::Integer;
                              left_tp::Union{Nothing, Tuple{Real, Real}}=nothing, right_tp::Union{Nothing, Tuple{Real, Real}}=nothing,
                              min_pos::Union{Nothing, Tuple{Real, Real}}=nothing, max_pos::Union{Nothing, Tuple{Real, Real}}=nothing,
                              turnpoint::Union{Real, Nothing}=nothing)
    𝑈 = x -> 𝐻₀(0.0, x, params)
    if left_tp === nothing # turning point intervals not provided, proceed to determine them
        left_tp, right_tp = turning_point_intervals(𝑈, min_pos, max_pos, turnpoint)
    end
    𝐸, 𝐸′, 𝐸″ = make_action_functions(𝑈, left_tp, right_tp)
    ClassicalHamiltonian(𝐻₀, 𝐻, 𝑈, Float64.(left_tp), Float64.(right_tp), 𝐸, 𝐸′, 𝐸″, params, s)
end

"""
Return the possible intervals of the turning points for motion in the potential 𝑈. A minimum and a maximum of the potential will be found
using the bracketing intervals `min_pos` and `max_pos`. If the heights of the "walls" of the potential are not equal, a `turnpoint` has to be provided.
For example, if the left wall is higher than the right one, the left turning point will be searched for in the interval (`turnpoint`, `x_min`).
"""
function turning_point_intervals(𝑈::Function, min_pos::Tuple{<:Real, <:Real}, max_pos::Tuple{<:Real, <:Real}, turnpoint::Union{Real, Nothing})
    # find position of the potential minimum
    result = Optim.optimize(x -> 𝑈(first(x)), min_pos[1], min_pos[2], Optim.Brent())
    x_min = Optim.minimizer(result)
    
    # find position of the potential maximum
    result = Optim.optimize(x -> -𝑈(first(x)), max_pos[1], max_pos[2], Optim.Brent())
    x_max = Optim.minimizer(result)
    E_max = -Optim.minimum(result)

    if x_max > x_min # if the located maximum is to the right of the minimum
        right_tp = (x_min, x_max)
        # if the `turnpoint` is not provided, calculate the left turning point assuming a symmetric well; 
        # otherwise, find the coordinate of the point giving `E_max` on the left wall
        x_max_left = turnpoint === nothing ? x_min - (x_max - x_min) :
                                             Roots.find_zero(x -> 𝑈(x) - E_max, (turnpoint, x_min), Roots.A42(), xrtol=1e-3)
        left_tp = (x_max_left, x_min)
    else # if the located maximum is to the left of the minimum
        left_tp = (x_max, x_min)
        x_max_right = turnpoint === nothing ? x_min + (x_min - x_max) :
                                              Roots.find_zero(x -> 𝑈(x) - E_max, (x_min, turnpoint), Roots.A42(), xrtol=1e-3)
        right_tp = (x_min, x_max_right)
    end
    left_tp, right_tp
end

"Construct and return the functions 𝐸(𝐼), 𝐸′(𝐼), and 𝐸″(𝐼)."
function make_action_functions(𝑈::Function, left_tp::Tuple{<:Real, <:Real}, right_tp::Tuple{<:Real, <:Real})
    n_E = 100 # number of energies (and actions) to save
    I = Vector{Float64}(undef, n_E) # for storing values of the action variable
    ΔU = 1e-4 * abs(𝑈(right_tp[2]) - 𝑈(right_tp[1])) # a slight shift for energies (see next line) to avoid problems during zero search
    E = range(𝑈(right_tp[1]) + ΔU, 𝑈(right_tp[2]) - ΔU, length=n_E) # energies inside the potential "well"

    for i in eachindex(E)
        x_min, x_max = turning_points(𝑈, E[i], left_tp, right_tp)
        # calculate ∫𝑝d𝑥 for a half-period; the second half is the same, hence no division by 2
        I[i] = quadgk(x -> 𝑝(𝑈, E[i], x), x_min, x_max, rtol=1e-4)[1] / π # `[1]` contains the integral, `[2]` contains error
    end
    
    𝐸  = Dierckx.Spline1D(I, E; k=2)      
    𝐸′ = x -> Dierckx.derivative(𝐸, x; nu=1)
    𝐸″ = x -> Dierckx.derivative(𝐸, x; nu=2)
    return 𝐸, 𝐸′, 𝐸″
end

"""
Return the turning points of motion with the given energy `E`. The initial guesses `a` and `b` should be given as
tuples representing the bracketing intervals.
"""
function turning_points(𝑈::Function, E::Real, a::Tuple{<:Real, <:Real}, b::Tuple{<:Real, <:Real})
    x_min = Roots.find_zero(x -> 𝑈(x) - E, a, atol=1e-4, Roots.A42())
    x_max = Roots.find_zero(x -> 𝑈(x) - E, b, atol=1e-4, Roots.A42())
    return x_min, x_max
end

"Momentum 𝑝(𝑥) = √[𝐸 - 𝑈(𝑥)] of a particle of energy `E`."
function 𝑝(𝑈::Function, E::Real, x::Real)
    p = E - 𝑈(x)
    p < 0 ? zero(p) : sqrt(p) # a safeguard for the case when `x` is slightly outside of the accessible region of oscillations
end

"""
Return action for the given energy `E` as the integral of momentum over a period of motion.
The turning points will be determined using the bracketing intervals `H.left_tp` and `H.right_tp`.
"""
function 𝐼(H::ClassicalHamiltonian, E::Real)
    x_min, x_max = turning_points(H.𝑈, E, H.left_tp, H.right_tp)
    # calculate ∫𝑝d𝑥 for a half-period; the second half is the same, hence no division by 2
    return quadgk(x -> 𝑝(H.𝑈, E, x), x_min, x_max, rtol=1e-4)[1] / π # `[1]` contains the integral, `[2]` contains error
end

"""
Return the action and mass at the working point. Also return the 𝑚th Fourier coefficient for every function in `perturbations`,
where the integer numbers 𝑚 are specified in `m`. `perturbations` are the spatial functions that couple the temporal perturbations;
Their signature is `f(p, x) = ...`.
"""
function compute_parameters(H::ClassicalHamiltonian, perturbations::Vector{Function}, m::Vector{<:Integer})
    ω = H.params[end]
    Ω = ω / H.s # our choice of the oscillation frequency (of the unperturbed system)
    Iₛ::Float64 = Roots.find_zero(x -> H.𝐸′(x) - Ω, 0.8last(Dierckx.get_knots(H.𝐸)), atol=1e-5) # find which 𝐼ₛ gives the frequency Ω
    E₀::Float64 = H.𝐸(Iₛ)     # energy of the system oscillating at the frequency Ω
    M::Float64 = 1 / H.𝐸″(Iₛ) # "mass" of the system oscillating at the frequency Ω
    # evolve the unperturbed system for one period 
    T = 2π / Ω
    tspan = (0.0, T)
    # initial conditions; they may be chosen arbitrary as long as the total energy equals `E₀`
    x₀ = H.right_tp[1]; p₀ = 𝑝(H.𝑈, E₀, x₀); # we choose position at the minimum and calculate the momentum
    H₀_problem = HamiltonianProblem(H.𝐻₀, p₀, x₀, tspan, H.params)
    dt = 2e-4
    sol = DiffEq.solve(H₀_problem, DiffEq.McAte3(); dt)

    # calculate the requested Fourier coefficient for every function in `perturbations`
    coeffs = Vector{ComplexF64}(undef, length(perturbations))
    V = Vector{Float64}(undef, length(sol.t)) # for storing perturbation evaluated in the solution points
    for (i, 𝑉) in enumerate(perturbations)
        V .= 𝑉.(sol[1, :], sol[2, :])
        coeffs[i] = fourier_coeff(V, m[i], dt, T)
    end
    return Iₛ, M, coeffs
end

"Calculate the `n`th complex Fourier coefficient of `f`. Simple trapezoid rule is used."
function fourier_coeff(f::AbstractVector, n::Int, dt::AbstractFloat, T::AbstractFloat)
    (sum(f[i] * cispi(-2n*(i-1)*dt/T) for i = eachindex(f)) - (f[1] + f[end])/2) * dt/T
end

"""
Compute evolutions (using the perturbed Hamiltonian) of 𝑝(𝑡) and 𝑥(𝑡) for the energy corresponding to `I_target`, and initial coordinate determined by `χ₀`, with |χ₀| ≤ 1.
For χ₀ > 0, the initial coordinate is 𝑥₀ + χ₀(𝑥ᵣ - 𝑥₀), and for χ₀ < 0, it is 𝑥₀ + χ₀(𝑥₀ - 𝑥ₗ),
where 𝑥₀ is the equilibrium point, while 𝑥ₗ and 𝑥ᵣ are left and right turning points for energy corresponding to `I_target`.
The pairs (𝑝(𝑡), 𝑥(𝑡)) are registered stroboscopically at the intervals of the period of external driving; `n_T` pairs are registered.

Then, transform the obtained (𝑝, 𝑥) pairs to (𝐼, ϑ) and return the results as a tuple of two vectors.
Transformation is performed as follows: for each pair (𝑝ᵢ, 𝑥ᵢ), the energy of the unperturbed motion is calculated as
𝐸ᵢ = 𝐻₀(𝑝ᵢ, 𝑥ᵢ), and the energy is then converted to action using the function 𝐼(𝐸).
To find the phase ϑᵢ, a period 𝑇ᵢ of unperturbed motion with energy 𝐸ᵢ is calculated, and the time moment 𝑡 corresponding to 
the pair (𝑝ᵢ, 𝑥ᵢ) is found. The phase is then given by ϑᵢ = 2π𝑡/𝑇ᵢ. Alternatively, a function converting the point (𝑝ᵢ, 𝑥ᵢ) to angle can be provided
as `point_to_angle(p, x, E, T) = ...`. This is useful if analytical solution of unperturbed motion is available.
Note that some energy 𝐸ⱼ may be such large (due to the perturbation) that the system is no longer confined to a single potential well. In that case,
no corresponding action 𝐼(𝐸ⱼ) exists. This will happen if `I_target` is too large. In that case, an info message will be printed,
and energies starting with 𝐸ⱼ will be ignored.
"""
function compute_IΘ(H::ClassicalHamiltonian, I_target::Real; χ₀::Real=0, n_T::Integer=100, point_to_angle::Union{Function, Nothing}=nothing)
    abs(χ₀) > 1 && begin @warn "|χ₀| ≤ 1 not satisfied. Setting χ₀ to 0."; χ₀ = 0 end
    
    ω = H.params[end]
    T_external = 2π / ω # period of the external driving
    tspan = (0.0, n_T * T_external)
    if χ₀ == 0
        q₀ = H.right_tp[1] # set iniital coordinate to the potential minimum (this position with positive momenutm defines the zero phase)
    elseif χ₀ > 0
        right_tp::Float64 = Roots.find_zero(x -> H.𝐻₀(0, x, H.params) - H.𝐸(I_target), H.right_tp) # right turning point for action `I_target`
        q₀ = H.right_tp[1] + χ₀ * (right_tp - H.right_tp[1])
    elseif χ₀ < 0
        left_tp::Float64 = Roots.find_zero(x -> H.𝐻₀(0, x, H.params) - H.𝐸(I_target), H.left_tp) # left turning point for action `I_target`
        q₀ = H.right_tp[1] + χ₀ * (H.right_tp[1] - left_tp) # note that `χ₀` is negative here
    end
    p₀::Float64 = 𝑝(H.𝑈, H.𝐸(I_target), q₀)

    H_problem = HamiltonianProblem(H.𝐻, p₀, q₀, tspan, H.params)
    sol = DiffEq.solve(H_problem, DiffEq.KahanLi8(); dt=2e-4, saveat=T_external)
    p = @view sol[1, :]
    x = @view sol[2, :]
    
    # Calculate the energies that the free system would possess if it was at `x` with momenta `p`
    E = Float64[]
    sizehint!(E, length(sol.t))
    for (pᵢ::Float64, xᵢ::Float64) in zip(p, x)
        Eᵢ::Float64 = H.𝐻₀(pᵢ, xᵢ, H.params)
        if Eᵢ < H.𝑈(H.left_tp[1])
            push!(E, Eᵢ)
        else # Interrupt if energy `Eᵢ` exceeds that of the barrier. All subsequent energies are of no interest then.
            @info "Perturbation deconfines the particle if it starts at action 𝐼 = $I_target."
            break
        end
    end

    I::Vector{Float64} = map(x -> 𝐼(H, x), E)
    Θ = similar(I)
    t::Float64 = 0 # initialise with a type to prevent `Any`

    if point_to_angle === nothing
        # for all the equations below, the initial position is chosen to be the potential minimum
        x₀ = H.right_tp[1]

        # find phases from the coordinates
        for i in eachindex(Θ)
            T_free::Float64 = 2π / H.𝐸′(I[i]) # period of the unperturbed motion at action `I[i]`
            tspan = (0.0, 1.02T_free) # take slightly more than `T_free`. Due to solver inaccuracies we might not get a full perdiod, and subsequent search will fail
            p₀ = 𝑝(H.𝑈, E[i], x₀)
            H₀_problem = HamiltonianProblem(H.𝐻₀, p₀, x₀, tspan, H.params)
            sln = DiffEq.solve(H₀_problem, DiffEq.McAte3(); dt=2e-4)

            # Find the time point when the equilibrium point x₀ (i.e. the potential minimum) is reached.
            # The coordinate will be greater than x₀ at times in (0; t_eq) and less than x₀ at times in (t_eq; T_free).
            t_eq::Float64 = Roots.find_zero(t -> sln(t)[2] - x₀, T_free/2)

            # If the coordinate `x[i]` is very close to potential minimum `x₀`, the momentum `p[i]` may lie just outside of the bracketing interval,
            # causing the root finding to fail. However, in that case `p[i]` is either very close to its maximum, meaning `t = 0`,
            # or is very close to the minimum, meaning `t = t_eq`. The two cases can be discerned by the sign of the momentum.
            if isapprox(x[i], x₀, atol=5e-3)
                t = p[i] > 0 ? 0.0 : t_eq
            else
                # use the sign of the coordinate to determine which part of the period the point (x[i]; p[i]) is in
                bracket = x[i] > x₀ ? (0.0, t_eq) : (t_eq, T_free)
                # Find the time corresponding to momentum `p[i]`:
                f = t -> sln(t)[1] - p[i] # construct the function whose root will be searched for
                # Check that `bracket` is indeed a bracketing interval. This might not be the case due to various inaccuracies.
                if prod(f.(bracket)) < 0
                    t = Roots.find_zero(f, bracket, Roots.A42(), xrtol=1e-3)
                else # otherwise, use the midpoint of the `bracket` as a starting point.
                    t = Roots.find_zero(f, (bracket[1]+bracket[2])/2) # Note that in this case the algorithm may occasionally converge to the zero in the wrong half of the period
                end
            end
            Θ[i] = 2π * t / T_free
        end
    else
        for i in eachindex(Θ)
            T_free::Float64 = 2π / H.𝐸′(I[i]) # period of the unperturbed motion at action `I[i]`
            Θ[i] = point_to_angle(p[i], x[i], E[i], T_free)
        end
    end
    return I, Θ
end

export ClassicalHamiltonian,
    make_action_functions,
    compute_parameters,
    compute_IΘ

end