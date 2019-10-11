| Testing  | Coverage | Documentation |
| :-----:  | :------: | :-----------: |
| [![Build Status](https://travis-ci.org/sisl/GaussianFilters.jl.svg?branch=master)](https://travis-ci.org/sisl/GaussianFilters.jl) | [![Coverage Status](https://coveralls.io/repos/github/sisl/GaussianFilters.jl/badge.svg?branch=master)](https://coveralls.io/github/sisl/GaussianFilters.jl?branch=master) |  [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://sisl.github.io/GaussianFilters.jl/latest) |

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

See documentation and examples for more details.

## Examples

Examples notebooks can be found in the `notebooks` folder:

[Kalman Filter Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/KF_2DMotionExample.ipynb)

[Extended Kalman Filter Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/EKF_SpinningSatelliteExample.ipynb)

[Unscented Kalman Filter Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/UKF_NonholonomicRobot.ipynb)

[GM-PHD Object Surveillance Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/GMPHD_SurveillanceExample.ipynb)

[GM-PHD Aircraft Carrier Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/GMPHD_AircraftCarrierExample.ipynb)
