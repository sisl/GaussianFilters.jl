__precompile__(true)
module GaussianFilters

using LinearAlgebra
using ForwardDiff
using Random
import Random: rand

# Kalman, Extended Kalman, Unscented Kalman Filters

export
	DynamicsModel,
	LinearDynamicsModel,
	NonlinearDynamicsModel,
	ObservationModel,
	LinearObservationModel,
	NonlinearObservationModel,
	predict,
	measure
include("models.jl")


export
	AbstractFilter,
	KalmanFilter,
	ExtendedKalmanFilter,
	UnscentedKalmanFilter,
	GaussianBelief

include("kf_classes.jl")

export 
	simulate_step,
	simulation,
	run_filter,
	unpack

include("simulate.jl")

export	
	update,
	update_with_info

include("kf.jl")
include("ekf.jl")

export
	unscented_transform,
	unscented_transform_inverse
include("ukf.jl")

# Gaussian Mixture PHD Filter

export
	Measurement,
	Dynamics,
	GaussianMixture,
	Spawn,
	PHDFilter
include("gmphd_classes.jl")

export
	update,
	predict,
	measure,
	prune
include("gmphd.jl")

export
	multiple_target_state_extraction
include("gmphd_extraction.jl")

# TODO: Add Labeled Multi-Bernoulli Filter

# Utilities
export
	belief_ellipse
include("utils.jl")

# POMDPs.jl integration stub. The actual implementation lives in
# ext/GaussianFiltersPOMDPsExt.jl and overrides this when POMDPs is
# loaded. The stub errors with a helpful message otherwise.
"""
    pomdps_updater(filter::AbstractFilter)

Wrap a GaussianFilters filter in a `POMDPs.Updater` so it can be passed
to POMDPs.jl simulators like `HistoryRecorder`. Requires POMDPs.jl to
be loaded; the implementation lives in a package extension.
"""
function pomdps_updater end

export pomdps_updater

end # module
