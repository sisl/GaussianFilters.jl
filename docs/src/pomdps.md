# POMDPs.jl Integration

```@meta
CurrentModule = GaussianFilters
```

GaussianFilters ships with a [POMDPs.jl](https://github.com/JuliaPOMDP/POMDPs.jl)
package extension that activates automatically when both packages are loaded.
The extension wires `AbstractFilter` into the POMDPs belief-updater interface
so that a Kalman, Extended Kalman, or Unscented Kalman filter can be used
directly as a `POMDPs.Updater` — including with simulators such as
`HistoryRecorder` and policy/planner code that expects the standard POMDPs
interface.

## Quick start

```julia
using GaussianFilters
using POMDPs
using POMDPTools

# 1. Build a filter as usual
dmodel = NonlinearDynamicsModel(f, W)
omodel = NonlinearObservationModel(h, V)
ekf    = ExtendedKalmanFilter(dmodel, omodel)

# 2. Wrap it as a POMDPs.Updater
updater = pomdps_updater(ekf)

# 3. Pass it to any POMDPs.jl simulator
hist = simulate(HistoryRecorder(max_steps=60), pomdp, policy, updater)
```

The wrapper is needed because POMDPs.jl simulators dispatch on the abstract
type `POMDPs.Updater`, and `AbstractFilter` cannot subtype `POMDPs.Updater`
directly without making POMDPs.jl a hard dependency.

## Direct dispatch (without the wrapper)

For simpler use cases that don't need the full `POMDPs.simulate` machinery,
the extension also overloads `POMDPs.update` and `POMDPs.initialize_belief`
directly on `AbstractFilter`:

```julia
b1 = POMDPs.update(ekf, b0, action, observation)
```

This is convenient when integrating with code that calls `POMDPs.update`
generically but does not require the `<:Updater` subtype constraint.

## Initial beliefs from distributions

The extension supports initializing a `GaussianBelief` from any multivariate
normal distribution. This is useful when the initial state of a `POMDP` is
given as an `MvNormal`:

```julia
using Distributions
prior = MvNormal([0.0, 0.0], [1.0 0.0; 0.0 0.5])
b0    = POMDPs.initialize_belief(updater, prior)
```

## API

```@docs
GaussianFilters.pomdps_updater
```

## Examples

Two example scripts live in the
[`examples/`](https://github.com/sisl/GaussianFilters.jl/tree/master/examples)
directory:

- [`pomdps_integration.jl`](https://github.com/sisl/GaussianFilters.jl/blob/master/examples/pomdps_integration.jl) —
  a minimal demonstration of using a `KalmanFilter` through `POMDPs.update`.

- [`pendulum_ekf_ilqr.jl`](https://github.com/sisl/GaussianFilters.jl/blob/master/examples/pendulum_ekf_ilqr.jl) —
  closed-loop stabilization of a noisy inverted pendulum observed only
  through its angular velocity. Defines a `PendulumPOMDP`, an
  `ILQRPolicy <: POMDPs.Policy` that runs iterative LQR on the belief
  mean (certainty equivalent control), wraps the EKF with `pomdps_updater`,
  and drives the whole closed loop through `POMDPs.simulate`.
