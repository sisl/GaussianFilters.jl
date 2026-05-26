# Examples

Self-contained example scripts demonstrating the filters in `GaussianFilters.jl`.

| Script | Filter | Description |
|---|---|---|
| `kf_2d_motion.jl` | KF | 2D point-mass with linear dynamics and measurement |
| `ekf_spinning_satellite.jl` | EKF | Nonlinear rigid-body rotation with saturated observations |
| `ukf_nonholonomic_robot.jl` | UKF | Differential-drive robot with range observation |
| `gmphd_surveillance.jl` | GM-PHD | Multi-target tracking in a surveillance region |
| `gmphd_aircraft_carrier.jl` | GM-PHD | Aircraft carrier scenario with birth/spawn models |
| `gmphd_tests.jl` | GM-PHD | Visualization helpers and pruning/merging demos |
| `pomdps_integration.jl` | KF | Uses a `KalmanFilter` as a `POMDPs.jl` belief updater |
| `pendulum_ekf_ilqr.jl` | EKF | Closed-loop stabilization of an inverted pendulum from partial (angle-only) observations. Defines a `PendulumPOMDP` and an `ILQRPolicy <: POMDPs.Policy`, wraps the EKF with `pomdps_updater`, and runs everything through `POMDPs.simulate`. Writes an animation to `examples/outputs/pendulum_ekf_ilqr.gif` (gitignored). |

## Running

The examples use a separate Julia project that depends on the local
`GaussianFilters` source via `[sources]` in `examples/Project.toml`.

```sh
julia --project=examples examples/kf_2d_motion.jl
```

Or interactively:

```julia
julia> ]
(@v1.11) pkg> activate examples
(examples) pkg> instantiate
julia> include("examples/kf_2d_motion.jl")
```
