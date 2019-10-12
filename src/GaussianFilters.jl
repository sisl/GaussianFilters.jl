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

export	update
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

end # module
