# Kalman-Class Filters

```@meta
CurrentModule = GaussianFilters
```

The Kalman, Extended Kalman, and Unscented Kalman filters are used to estimate state using unimodal multivariate Gaussian distributions. A state estimate in this framework is defined as a `GaussianBelief` consisting of a mean and covariance.

```@docs
GaussianFilters.GaussianBelief
```

## Building a Filter

In general, Kalman-class filters can be built with either linear or non-linear dynamics and measurement models. Linear models should be defined with appropriately matrices. Non-linear models should be defined using an appropriate function of two variables, state and action. Both models should be defined with symmetric noise covariance matrices.

NOTE: There is no need to define Jacobians for non-linear models, since this package uses automatic forward differentiation to compute Jacobians in real time. Just make sure the models are forward differentiable in all possible belief locations.

```@docs
GaussianFilters.LinearDynamicsModel
GaussianFilters.LinearObservationModel
GaussianFilters.NonlinearDynamicsModel
GaussianFilters.NonlinearObservationModel
```

Use the filter constructors with the appropriately typed models to build a filter. It is recommended to always construct a Kalman filter type when both dynamics and observation models are linear.

```@docs
GaussianFilters.KalmanFilter
GaussianFilters.ExtendedKalmanFilter
GaussianFilters.UnscentedKalmanFilter
```


## Simulating Data

Given a filter, an initial belief, and an action sequence, you can either simulate state and measurement data all at once with `simulation` or one step at a time with `simulate_step`

```@docs
GaussianFilters.simulation
GaussianFilters.simulate_step
```

In addition, the dynamics and observation models can be queried on a single state control input using the `predict` and `measure` methods respectively.

```@docs 
    predict(::LinearDynamicsModel, ::AbstractVector{<:Number}, ::AbstractVector{<:Number})
    predict(::NonlinearDynamicsModel, ::AbstractVector{<:Number}, ::AbstractVector{<:Number})
    measure(::LinearObservationModel, ::AbstractVector{<:Number}, ::AbstractVector{<:Number})
    measure(::NonlinearObservationModel, ::AbstractVector{<:Number}, ::AbstractVector{<:Number})
```

## Running a Filter

You can run a filter on a sequential measurement data using the `run_filter` function.

```@docs
GaussianFilters.run_filter
```

Alternatively, you can make step-wise belief updates using the `update` function, which consists of a two-step process to a) `predict` the next state given a known action and b) make measurement-based belief updates with `measure`.

```@docs
GaussianFilters.update
GaussianFilters.predict
GaussianFilters.measure
```


## Utilities

The output of running a filter is a `GaussianBelief` vector, which can be condensed into appropriate tensors with `unpack`.

```@docs
GaussianFilters.unpack
```

`belief_ellipse` can be used to convert a 2-D Gaussian belief into points along a confidence interval ellipse for plotting.

```@docs
GaussianFilters.belief_ellipse
```

## Examples

Full implementation examples can be found in the `notebooks` folder of the repo:

[Kalman Filter Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/KF_2DMotionExample.ipynb)

[Extended Kalman Filter Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/EKF_SpinningSatelliteExample.ipynb)

[Unscented Kalman Filter Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/UKF_NonholonomicRobot.ipynb)
