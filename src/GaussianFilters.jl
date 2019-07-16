__precompile__(true)
module GaussianFilters


using LinearAlgebra
using Random
import Random: rand
import Base: step

# Kalman, Extended Kalman, Unscented Kalman Filters

export
	DynamicsModel,
	LinearDynamicsModeal,
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
	measure
include("kf.jl")

include("ekf.jl")

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

# Utilities

export
	beliefEllipse
include("utils.jl")

end # module
