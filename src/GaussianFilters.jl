__precompile__(true)
module GaussianFilters

using LinearAlgebra
using ForwardDiff
using Random
import Random: rand
import Base: step

# Kalman, Extended Kalman, Unscented Kalman Filters

export
	DynamicsModel,
	LinearDynamicsModel,
	NonlinearDynamicsModel,
	ObservationModel,
	LinearObservationModel,
	NonlinearObservationModel,
	AbstractFilter,
	KalmanFilter,
	ExtendedKalmanFilter,
	UnscentedKalmanFilter,
	GaussianBelief
include("kf_classes.jl")

export
	update,
	predict,
	measure,
	simulation,
	simulate_step,
	run_filter,
	unpack
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
	step,
	step_prune
include("phd_step.jl")

export
	prune
include("prune.jl")

export
	multiple_target_state_extraction
include("extraction.jl")

# TODO: Add Labeled Multi-Bernoulli Filter

# Utilities

export
	belief_ellipse
include("utils.jl")

end # module
