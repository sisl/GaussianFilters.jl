# Stabilizing a noisy inverted pendulum from partial observations
# ================================================================
#
# Demonstrates the GaussianFilters <-> POMDPs.jl integration in a closed
# control loop. At each step we
#   1. observe a noisy angle (the angular velocity is hidden),
#   2. update the EKF belief via POMDPs.update,
#   3. plan a torque sequence with iLQR using belief.μ as the starting
#      state (certainty equivalent control — iLQR ignores belief.Σ),
#   4. apply the first planned action.
#
# The dynamics are pulled from POMDPModels.InvertedPendulum so we are
# integrating against the real ecosystem rather than a private copy.

using GaussianFilters
using POMDPs
using POMDPModels: InvertedPendulum
using LinearAlgebra
using Random
using ForwardDiff
using Distributions: MvNormal

# ----- dynamics & observation models -----
#
# We pull parameters from POMDPModels.InvertedPendulum and re-implement
# the Euler step generically (the POMDPModels.euler signature pins state
# to Tuple{Float64,Float64}, which blocks ForwardDiff dual numbers).

const ip = InvertedPendulum()
const σ_proc = 0.01    # process noise std per state component
const σ_obs  = 0.05    # angle observation noise std (rad)

function dwdt(th, w, u, m=ip)
    num = m.g*sin(th) - m.alpha*m.m*m.l*(w^2)*sin(2*th)*0.5 - m.alpha*cos(th)*u
    den = (4/3)*m.l - m.alpha*m.l*(cos(th)^2)
    return num/den
end

function pendulum_step(x::AbstractVector, u::AbstractVector)
    th, w, ut = x[1], x[2], u[1]
    alph = dwdt(th, w, ut)
    return [th + w*ip.dt + 0.5*alph*ip.dt^2,
            w + alph*ip.dt]
end

# Observation: noisy measurement of the angle only.
observe_angle(x::AbstractVector, u::AbstractVector) = [x[1]]

W = σ_proc^2 * Matrix{Float64}(I, 2, 2)
V = (σ_obs^2)  * Matrix{Float64}(I, 1, 1)

dmodel = NonlinearDynamicsModel(pendulum_step, W)
omodel = NonlinearObservationModel(observe_angle, V)
ekf = ExtendedKalmanFilter(dmodel, omodel)

# ----- minimal iLQR -----
#
# Standard LQR backward pass on the linearization of pendulum_step around
# the current nominal trajectory, with a quadratic stage cost
#
#   ℓ(x,u) = xᵀ Q x + uᵀ R u
#
# and quadratic terminal cost xᵀ Qf x. We linearize via ForwardDiff at
# every step of the forward roll, do one backward LQR sweep to compute
# feedback gains, and return the resulting control sequence. No line
# search, no regularization — sufficient to stabilize from a small
# perturbation (θ₀ ~ 0.3) but not to swing up from θ = π.

const Q  = Diagonal([10.0, 1.0])        # stage cost on (θ, ω)
const Qf = Diagonal([100.0, 10.0])      # terminal cost
const R  = Diagonal([0.01])             # control effort

function rollout(x0::Vector{Float64}, us::Vector{Vector{Float64}})
    H = length(us)
    xs = Vector{Vector{Float64}}(undef, H + 1)
    xs[1] = x0
    for t in 1:H
        xs[t+1] = pendulum_step(xs[t], us[t])
    end
    return xs
end

function ilqr_step(x0::Vector{Float64}; H::Int = 20)
    # Initialize nominal control sequence at zero, roll out, then take
    # one LQR backward pass on the linearization of that trajectory.
    us = [zeros(1) for _ in 1:H]
    xs = rollout(x0, us)

    # Terminal value
    Vx  = 2 * Qf * xs[end]
    Vxx = 2 * Qf

    K = Vector{Matrix{Float64}}(undef, H)
    k = Vector{Vector{Float64}}(undef, H)

    for t in H:-1:1
        x, u = xs[t], us[t]
        # Jacobians of dynamics at (x, u)
        fx = ForwardDiff.jacobian(z -> pendulum_step(z, u), x)
        fu = ForwardDiff.jacobian(z -> pendulum_step(x, z), u)

        # Stage cost derivatives
        lx  = 2 * Q * x
        lu  = 2 * R * u
        lxx = 2 * Q
        luu = 2 * R

        # Q-function expansion
        Qx  = lx  + fx' * Vx
        Qu  = lu  + fu' * Vx
        Qxx = lxx + fx' * Vxx * fx
        Quu = luu + fu' * Vxx * fu
        Qux = fu' * Vxx * fx

        # Feedback gains
        Quu_inv = inv(Quu)
        k[t] = -Quu_inv * Qu
        K[t] = -Quu_inv * Qux

        # Value-function update
        Vx  = Qx + K[t]' * Quu * k[t] + K[t]' * Qu + Qux' * k[t]
        Vxx = Qxx + K[t]' * Quu * K[t] + K[t]' * Qux + Qux' * K[t]
    end

    # Forward roll of the new control law (no line search).
    new_us = Vector{Vector{Float64}}(undef, H)
    new_xs = Vector{Vector{Float64}}(undef, H + 1)
    new_xs[1] = x0
    for t in 1:H
        new_us[t] = us[t] + k[t] + K[t] * (new_xs[t] - xs[t])
        new_xs[t+1] = pendulum_step(new_xs[t], new_us[t])
    end
    return new_us
end

# ----- closed-loop simulation -----

rng = MersenneTwister(0)

# True initial state: a small angle perturbation from upright.
x_true = [0.3, 0.0]

# Initial belief: slightly uncertain prior centered near the true state.
prior = MvNormal([0.3, 0.0], Matrix(Diagonal([0.05, 0.1])))
belief = POMDPs.initialize_belief(ekf, prior)

T = 60                       # control steps (60 * dt = 6 s)
u_prev = [0.0]
angle_history = Float64[x_true[1]]
mean_history = [copy(belief.μ)]

for t in 1:T
    global x_true, belief, u_prev

    # Plan a torque from the certainty-equivalent state.
    us = ilqr_step(belief.μ)
    u  = us[1]

    # Step the true system (with process noise).
    x_true = pendulum_step(x_true, u) .+ σ_proc .* randn(rng, 2)

    # Observe (with measurement noise).
    y = observe_angle(x_true, u) .+ σ_obs .* randn(rng, 1)

    # Belief update via the POMDPs.jl interface.
    belief = POMDPs.update(ekf, belief, u, y)

    push!(angle_history, x_true[1])
    push!(mean_history, copy(belief.μ))
    u_prev = u
end

θf = angle_history[end]
ωf = mean_history[end][2]
println("Final true angle:        ", round(θf, digits=4), " rad")
println("Final belief mean (θ,ω): (",
        round(mean_history[end][1], digits=4), ", ",
        round(ωf, digits=4), ")")
println("Max |angle| over run:    ", round(maximum(abs, angle_history), digits=4), " rad")
