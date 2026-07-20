# Stabilizing a noisy inverted pendulum with EKF beliefs + iLQR control
# =====================================================================
#
# This example defines a partial-observation `PendulumPOMDP`, uses an
# `ExtendedKalmanFilter` wrapped as a `POMDPs.Updater`, and an
# `ILQRPolicy` that runs iterative LQR on the belief mean (certainty
# equivalence — iLQR ignores belief covariance). The closed loop is
# driven by `POMDPs.simulate`, not a hand-rolled control loop.
#
# The observation is the angular velocity ω only (gyroscope-style);
# the angle θ is hidden and must be inferred. This is the more
# pedagogically interesting filtering problem: the controller acts on
# `belief.μ[1]` for the angle even though the angle is never directly
# observed. The animation shows σθ shrinking over time (the filter
# converging) while σω stays small throughout (direct observation).

using GaussianFilters
using POMDPs
using POMDPTools: HistoryRecorder, eachstep
using POMDPModels: InvertedPendulum
using LinearAlgebra
using Random
using ForwardDiff
using StaticArrays
import Distributions
using Distributions: MvNormal

# Plots is needed for the animation block at the end. Loading it up
# front so that macros (@animate, gif) are available at parse time.
ENV["GKSwstype"] = "100"
using Plots

# ---------------------------------------------------------------------
# Dynamics
# ---------------------------------------------------------------------
# Parameters are pulled from POMDPModels.InvertedPendulum so we are not
# inventing physical constants. We re-express the Euler step generically
# because POMDPModels.euler is signature-locked to Tuple{Float64,Float64}
# and blocks ForwardDiff dual numbers (which are used both inside the
# EKF Jacobian and inside the iLQR linearization).

const IP = InvertedPendulum()

function dwdt(th, w, u, m = IP)
    num = m.g*sin(th) - m.alpha*m.m*m.l*(w^2)*sin(2*th)*0.5 - m.alpha*cos(th)*u
    den = (4/3)*m.l - m.alpha*m.l*(cos(th)^2)
    return num / den
end

# Deterministic Euler step. Returns an SVector to exercise our
# StaticArrays compatibility end-to-end.
function pendulum_step(x::AbstractVector, u::AbstractVector)
    th, w, ut = x[1], x[2], u[1]
    alph = dwdt(th, w, ut)
    return SVector{2}(th + w*IP.dt + 0.5*alph*IP.dt^2,
                      w + alph*IP.dt)
end

# Observe the angular velocity only (gyroscope-like measurement). The
# angle itself must be inferred by the EKF — there is no direct
# observation of absolute position. This is the more pedagogically
# interesting variant: the controller acts on belief.μ for θ, which
# is reconstructed entirely from successive ω observations.
observe_omega(x::AbstractVector, u::AbstractVector) = SVector{1}(x[2])

# ---------------------------------------------------------------------
# POMDP definition
# ---------------------------------------------------------------------

const σ_proc = 0.01    # process noise std (per state component)
const σ_obs  = 0.1     # angular velocity observation noise std (rad/s)

struct PendulumPOMDP <: POMDP{SVector{2,Float64}, SVector{1,Float64}, SVector{1,Float64}}
    σ_proc::Float64
    σ_obs::Float64
    dt::Float64
    θ_max::Float64       # episode terminates if |θ| exceeds this
    γ::Float64
end

PendulumPOMDP() = PendulumPOMDP(σ_proc, σ_obs, IP.dt, π, 0.99)

# Initial-state distribution: a Gaussian over (θ, ω). HistoryRecorder
# samples the true initial state from this and `initialize_belief` reads
# the same distribution's mean/cov (via our extension) to seed the EKF.
# We wrap MvNormal so `rand` returns an SVector instead of a Vector.
struct SVecMvNormal
    mvn::MvNormal
end
Base.rand(rng::AbstractRNG, d::SVecMvNormal) = SVector{2}(rand(rng, d.mvn))
Distributions.mean(d::SVecMvNormal) = mean(d.mvn)
Distributions.cov(d::SVecMvNormal)  = cov(d.mvn)
# Forward initialize_belief to the underlying MvNormal so our extension picks it up.
POMDPs.initialize_belief(u::POMDPs.Updater, d::SVecMvNormal) = POMDPs.initialize_belief(u, d.mvn)

POMDPs.initialstate(::PendulumPOMDP) =
    SVecMvNormal(MvNormal([0.0, 0.0], Matrix(Diagonal([0.5^2, 0.1^2]))))

POMDPs.isterminal(p::PendulumPOMDP, s::AbstractVector) = abs(s[1]) > p.θ_max
POMDPs.discount(p::PendulumPOMDP) = p.γ

# `gen` is the generative interface: (next_state, observation, reward).
# Relaxing the state argument to AbstractVector (rather than SVector{2})
# accommodates POMDPs' internal handling.
function POMDPs.gen(p::PendulumPOMDP, s::AbstractVector, a::AbstractVector, rng)
    sp = pendulum_step(s, a) .+ p.σ_proc .* SVector{2}(randn(rng), randn(rng))
    o  = observe_omega(sp, a) .+ p.σ_obs .* SVector{1}(randn(rng))
    r  = -(s[1]^2 + 0.1*s[2]^2 + 0.01*a[1]^2)
    return (sp = sp, o = o, r = r)
end

# ---------------------------------------------------------------------
# ILQRPolicy
# ---------------------------------------------------------------------
# Holds only iLQR hyperparameters + warm-start solver state. No belief,
# no EKF — that's the Updater's job.

mutable struct ILQRPolicy{TQ, TR} <: POMDPs.Policy
    Q::TQ
    Qf::TQ
    R::TR
    H::Int                                  # planning horizon
    max_iters::Int
    u_max::Float64                          # actuator saturation
    us_warm::Vector{SVector{1,Float64}}     # solver warm start
end

function ILQRPolicy(;
        Q  = Diagonal(SVector{2}(10.0, 1.0)),
        Qf = Diagonal(SVector{2}(50.0, 5.0)),
        R  = Diagonal(SVector{1}(0.05)),
        H::Int = 40,
        max_iters::Int = 8,
        u_max::Float64 = 100.0)
    return ILQRPolicy(Q, Qf, R, H, max_iters, u_max, SVector{1,Float64}[])
end

# ---------------------------------------------------------------------
# iLQR free functions (take the policy as input)
# ---------------------------------------------------------------------

function rollout(x0, us, p::ILQRPolicy)
    H = length(us)
    xs = Vector{SVector{2,Float64}}(undef, H + 1)
    xs[1] = x0
    for t in 1:H
        xs[t+1] = pendulum_step(xs[t], us[t])
    end
    return xs
end

function trajectory_cost(xs, us, p::ILQRPolicy)
    c = xs[end]' * p.Qf * xs[end]
    for t in eachindex(us)
        c += xs[t]' * p.Q * xs[t] + us[t]' * p.R * us[t]
    end
    return c
end

function backward_pass(xs, us, p::ILQRPolicy, μ_reg)
    H = length(us)
    Vx  = 2 * p.Qf * xs[end]
    Vxx = 2 * Matrix(p.Qf)
    K = Vector{Matrix{Float64}}(undef, H)
    k = Vector{Vector{Float64}}(undef, H)
    for t in H:-1:1
        x, u = xs[t], us[t]
        fx = ForwardDiff.jacobian(z -> pendulum_step(z, u), Vector(x))
        fu = ForwardDiff.jacobian(z -> pendulum_step(x, z), Vector(u))
        Qx  = 2 * Vector(p.Q  * x) + fx' * Vx
        Qu  = 2 * Vector(p.R  * u) + fu' * Vx
        Qxx = 2 * Matrix(p.Q)      + fx' * Vxx * fx
        Quu = 2 * Matrix(p.R)      + fu' * Vxx * fu
        Qux = fu' * Vxx * fx
        Quu_inv = inv(Quu + μ_reg * I)
        k[t] = -Quu_inv * Qu
        K[t] = -Quu_inv * Qux
        Vx  = Qx + K[t]' * Quu * k[t] + K[t]' * Qu + Qux' * k[t]
        Vxx = Qxx + K[t]' * Quu * K[t] + K[t]' * Qux + Qux' * K[t]
    end
    return K, k
end

function forward_pass(xs, us, K, k, α, p::ILQRPolicy)
    H = length(us)
    new_xs = Vector{SVector{2,Float64}}(undef, H + 1)
    new_us = Vector{SVector{1,Float64}}(undef, H)
    new_xs[1] = xs[1]
    for t in 1:H
        u_lin = us[t] + α * SVector{1}(k[t][1]) + SVector{1}(K[t][1,1]*(new_xs[t][1]-xs[t][1]) +
                                                              K[t][1,2]*(new_xs[t][2]-xs[t][2]))
        new_us[t] = SVector{1}(clamp(u_lin[1], -p.u_max, p.u_max))
        new_xs[t+1] = pendulum_step(new_xs[t], new_us[t])
    end
    return new_xs, new_us
end

function ilqr(x0::SVector{2,Float64}, p::ILQRPolicy)
    us = isempty(p.us_warm) ? [SVector{1}(0.0) for _ in 1:p.H] : copy(p.us_warm)
    if length(us) != p.H
        us = [SVector{1}(0.0) for _ in 1:p.H]
    end
    xs = rollout(x0, us, p)
    J  = trajectory_cost(xs, us, p)
    μ_reg = 1e-3
    for _ in 1:p.max_iters
        K, k = backward_pass(xs, us, p, μ_reg)
        accepted = false
        α = 1.0
        for _ in 1:10
            new_xs, new_us = forward_pass(xs, us, K, k, α, p)
            J_new = trajectory_cost(new_xs, new_us, p)
            if J_new < J - 1e-4
                xs, us = new_xs, new_us
                J = J_new
                accepted = true
                μ_reg = max(μ_reg / 2, 1e-6)
                break
            end
            α *= 0.5
        end
        accepted || (μ_reg *= 4)
        μ_reg > 1e6 && break
    end
    return us
end

# ---------------------------------------------------------------------
# THE POMDPs.Policy interface — a one-liner.
# ---------------------------------------------------------------------
function POMDPs.action(policy::ILQRPolicy, belief::GaussianBelief)
    us = ilqr(SVector{2}(belief.μ[1], belief.μ[2]), policy)
    # Shift solution by one step so the next call warm-starts.
    policy.us_warm = vcat(us[2:end], [SVector{1}(0.0)])
    return us[1]
end

# ---------------------------------------------------------------------
# Build the components
# ---------------------------------------------------------------------

const W = (σ_proc^2) * Matrix{Float64}(I, 2, 2)
const V = (σ_obs^2)  * Matrix{Float64}(I, 1, 1)

dmodel = NonlinearDynamicsModel(pendulum_step, W)
omodel = NonlinearObservationModel(observe_omega, V)
ekf    = ExtendedKalmanFilter(dmodel, omodel)
updater = pomdps_updater(ekf)          # wraps the EKF as a POMDPs.Updater
pomdp   = PendulumPOMDP()
policy  = ILQRPolicy()

# ---------------------------------------------------------------------
# Run it via POMDPs.simulate
# ---------------------------------------------------------------------

rng = MersenneTwister(0)
hr  = HistoryRecorder(rng = rng, max_steps = 60)
hist = simulate(hr, pomdp, policy, updater)

# Pull trajectory out of the history for reporting (and downstream
# plotting, etc).
angle_history  = [step.s[1] for step in eachstep(hist)]
push!(angle_history, hist[end].sp[1])              # final state
omega_history  = [step.s[2] for step in eachstep(hist)]
push!(omega_history, hist[end].sp[2])
belief_history = [step.b for step in eachstep(hist)]
push!(belief_history, hist[end].bp)                # final belief
action_history = [step.a[1] for step in eachstep(hist)]

θ_final = angle_history[end]
θ_max   = maximum(abs, angle_history)
steady = angle_history[(2*length(angle_history))÷3 : end]
θ_ss_mean = sum(steady) / length(steady)
θ_ss_max  = maximum(abs, steady)

println("Initial angle:                          ", round(angle_history[1], digits=4), " rad")
println("Final true angle:                       ", round(θ_final, digits=4), " rad")
println("Max |angle| during run:                 ", round(θ_max, digits=4), " rad")
println("Steady-state mean |angle| (last third): ", round(abs(θ_ss_mean), digits=4), " rad")
println("Steady-state max  |angle| (last third): ", round(θ_ss_max, digits=4), " rad")
println()
println("Final belief mean (θ, ω):               (",
        round(belief_history[end].μ[1], digits=4), ", ",
        round(belief_history[end].μ[2], digits=4), ")")

# ---------------------------------------------------------------------
# Animation (written to examples/outputs/, which is gitignored)
# ---------------------------------------------------------------------
# Set GKSwstype=100 so Plots runs headless. Set GENERATE_GIF=false in
# the environment to skip rendering.

if get(ENV, "GENERATE_GIF", "true") != "false"
    out_dir = joinpath(@__DIR__, "outputs")
    mkpath(out_dir)
    gif_path = joinpath(out_dir, "pendulum_ekf_ilqr.gif")

    N = length(angle_history)
    # Pre-compute belief mean and ±2σ envelopes over the full run.
    μθ = [b.μ[1] for b in belief_history]
    μω = [b.μ[2] for b in belief_history]
    σθ = [2*sqrt(b.Σ[1,1]) for b in belief_history]
    σω = [2*sqrt(b.Σ[2,2]) for b in belief_history]

    anim = @animate for i in 1:N
        θ_true = angle_history[i]
        θ_blf  = belief_history[i].μ[1]

        # --- left: pendulum rod animation ---
        L = 1.0
        x_tip_true, y_tip_true = L*sin(θ_true), L*cos(θ_true)
        x_tip_blf,  y_tip_blf  = L*sin(θ_blf),  L*cos(θ_blf)

        p1 = plot([0.0, x_tip_true], [0.0, y_tip_true],
                  lw=4, color=:steelblue, label="true",
                  xlim=(-1.2, 1.2), ylim=(-0.4, 1.2),
                  aspect_ratio=:equal, framestyle=:box,
                  title="Pendulum (step $(i-1)/$(N-1))")
        plot!(p1, [0.0, x_tip_blf], [0.0, y_tip_blf],
              lw=2, color=:orange, linestyle=:dash, label="belief")
        scatter!(p1, [0.0], [0.0], color=:black, ms=4, label=false)
        scatter!(p1, [x_tip_true], [y_tip_true], color=:steelblue, ms=6, label=false)

        # --- top right: θ over time with belief ±2σ ribbon ---
        ts = 0:(i-1)
        p2 = plot(ts, μθ[1:i], ribbon=σθ[1:i],
                  color=:orange, lw=2, fillalpha=0.25,
                  label="belief μ ± 2σ",
                  xlim=(0, N-1), ylim=(-1.0, 1.0),
                  ylabel="θ (rad)", legend=:topright,
                  guidefonthalign=:right, ymirror=false,
                  left_margin=8Plots.mm)
        plot!(p2, ts, angle_history[1:i],
              color=:steelblue, lw=2, label="θ true")
        hline!(p2, [0.0], color=:gray, linestyle=:dot, label=false)

        # --- bottom right: ω over time with belief ±2σ ribbon ---
        p3 = plot(ts, μω[1:i], ribbon=σω[1:i],
                  color=:orange, lw=2, fillalpha=0.25,
                  label="belief μ ± 2σ",
                  xlim=(0, N-1), ylim=(-3.0, 3.0),
                  xlabel="step", ylabel="ω (rad/s)", legend=:topright,
                  guidefonthalign=:right,
                  left_margin=8Plots.mm, bottom_margin=8Plots.mm)
        plot!(p3, ts, omega_history[1:i],
              color=:steelblue, lw=2, label="ω true (observed)")
        hline!(p3, [0.0], color=:gray, linestyle=:dot, label=false)

        right = plot(p2, p3, layout=(2, 1))
        plot(p1, right, layout=(1, 2), size=(1100, 550),
             bottom_margin=5Plots.mm)
    end
    gif(anim, gif_path, fps=10)
    println()
    println("Animation saved to: $gif_path")
end
