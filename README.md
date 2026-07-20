[![CI](https://github.com/sisl/GaussianFilters.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/sisl/GaussianFilters.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/sisl/GaussianFilters.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/sisl/GaussianFilters.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://sisl.github.io/GaussianFilters.jl/latest)

# GaussianFilters.jl

GaussianFilters implements methods to define and run **Kalman**, **Extended Kalman**, **Unscented Kalman**, and **Gaussian-Mixture Probability Hypothesis Density** Filters on simulated data. It also implements simulation functions for the Kalman-class filters.

## Documentation

The documentation for the package can be found here: <https://sisl.github.io/GaussianFilters.jl/latest>

## Installation

GaussianFilters can be installed by running:

```julia
using Pkg
Pkg.add("GaussianFilters")
```

## Basic Usage

Basic usage follows along defining appropriate models, constructing an appropriate filter, and running the filter with known actions on some measurement data.

```julia
using GaussianFilters, LinearAlgebra

# dynamics model
A = [1 0.1; 0 1]
B = [0; 1]
W = [0.5 0; 0 0.5]
dmodel = LinearDynamicsModel(A, B, W)

# measurement model
measure(x, u) = LinearAlgebra.norm(x, 2)
V = [0.01]
omodel = NonlinearObservationModel(measure, V)

# filtering given some action and measurement
ukf = UnscentedKalmanFilter(dmodel, omodel)

b0 = GaussianBelief([0, 0], [1 0; 0 1])
b1 = update(ukf, b0, action, measurement)
```

Covariances may be passed as any `AbstractMatrix`, or in identity form
(`I(n)`, `s*I(n)`). `GaussianBelief` also accepts a bare `UniformScaling`,
inferring the size from the mean vector:

```julia
b0 = GaussianBelief([0.0, 0.0], 2.0*I)   # Σ becomes 2*I(2)
```

See the `examples/` directory for full demonstrations.

## Examples

Each example script lives in `examples/` with its own `Project.toml`.
Run any of them with:

```sh
julia --project=examples examples/kf_2d_motion.jl
```

| Script | Filter | Description |
|---|---|---|
| [`kf_2d_motion.jl`](examples/kf_2d_motion.jl) | KF | 2D point-mass with linear dynamics |
| [`ekf_spinning_satellite.jl`](examples/ekf_spinning_satellite.jl) | EKF | Rigid-body rotation with saturated observations |
| [`ukf_nonholonomic_robot.jl`](examples/ukf_nonholonomic_robot.jl) | UKF | Differential-drive robot with range observation |
| [`gmphd_surveillance.jl`](examples/gmphd_surveillance.jl) | GM-PHD | Multi-target tracking |
| [`gmphd_aircraft_carrier.jl`](examples/gmphd_aircraft_carrier.jl) | GM-PHD | Aircraft carrier scenario |
| [`gmphd_tests.jl`](examples/gmphd_tests.jl) | GM-PHD | Visualization helpers |
