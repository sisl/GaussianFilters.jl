# Gaussian Mixture Probability Hypothesis Density Filter

```@meta
CurrentModule = GaussianFilters
```

A Gaussian Mixture Probability Hypothesis Density (GM-PHD) Filter handles multi-object tracking in a low signal-to-noise environment by recursively propagating an unnormalized [Gaussian Mixture Model](https://en.wikipedia.org/wiki/Mixture_model#Multivariate_Gaussian_mixture_model), which when integrated over a volume corresponds with the expected number of objects in that volume.

```@docs
GaussianFilters.GaussianMixture
```

## Building a GM-PHD Filter

Building a GM-PHD Filter requires defining a birth intensity model of type `GaussianMixture` for where objects can be globally born, a spawn intensity model of type `Spawn` for how new objects can spawn off old objects, a Vector of possible linear dynamics models of type `Dynamics`, a linear measurement model of type `Measurement`, a survival probability, a detection probability, and a clutter modeling function.

```@docs
GaussianFilters.Spawn
GaussianFilters.Dynamics
GaussianFilters.Measurement
```

Once all of these are defined, a GM-PHD filter can be constructed with `PHDFilter`.

```@docs
GaussianFilters.PHDFilter
```

## Running a GM-PHD Filter

Currently, the GM-PHD Filter must be run step-by-step. Similar to the Kalman-Class filters, this can be done with a single call to `update`, which wraps functions to `predict` the next state, perform a measurement update with `measure`, and `prune` the resulting mixture model of low-probability and sufficiently close mixtures.

```@docs
update(::PHDFilter, ::GaussianMixture, ::Vector{<:AbstractVector{<:Number}}, ::Real, ::Real, ::Integer)
predict(::PHDFilter, ::GaussianMixture)
measure(::PHDFilter, ::GaussianMixture, ::Vector{<:AbstractVector{<:Real}})
GaussianFilters.prune
```

The GM-PHD update internally evaluates a multivariate normal density:

```@docs
GaussianFilters.MvNormalPDF
```

Target locations can be extracted from a `GaussianMixture` state using `multiple_target_state_extraction`.

```@docs
GaussianFilters.multiple_target_state_extraction
```

## Examples

Full implementation examples can be found in the [`examples/`](https://github.com/sisl/GaussianFilters.jl/tree/master/examples) directory of the repo:

- [`gmphd_surveillance.jl`](https://github.com/sisl/GaussianFilters.jl/blob/master/examples/gmphd_surveillance.jl) — Object Surveillance
- [`gmphd_aircraft_carrier.jl`](https://github.com/sisl/GaussianFilters.jl/blob/master/examples/gmphd_aircraft_carrier.jl) — Aircraft Carrier scenario
- [`gmphd_tests.jl`](https://github.com/sisl/GaussianFilters.jl/blob/master/examples/gmphd_tests.jl) — Visualization helpers
