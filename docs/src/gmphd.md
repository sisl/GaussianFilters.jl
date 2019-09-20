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

Currently, the GM-PHD Filter must be run step-by-step. Prediction and measurement steps are combined in `step`, while resulting low-probability mixture components are pruned with `prune`. These functions are combined in `step_prune`.

```@docs
GaussianFilters.step
GaussianFilters.prune
GaussianFilters.step_prune
```
Target locations can be extracted from a `GaussianMixture` state using `multiple_target_state_extraction`.

```@docs
GaussianFilters.multiple_target_state_extraction
```

## Examples

Full implementation examples can be found in the `notebooks` folder of the repo:

[GM-PHD Object Surveillance Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/GMPHD_SurveillanceExample.ipynb)

[GM-PHD Aircraft Carrier Example](https://github.com/sisl/GaussianFilters.jl/blob/master/notebooks/GMPHD_AircraftCarrierExample.ipynb)
