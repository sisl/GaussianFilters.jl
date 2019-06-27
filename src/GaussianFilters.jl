__precompile__(true)
module GaussianFilters


using LinearAlgebra
using Random
import Random: rand
import Base: step

# Kalman, Extended Kalman, Unscented Kalman Filters

# TODO #

# Gaussian Mixture PHD Filter

export
	Measurement,
	Dynamics,
	GaussianMixture,
	Spawn,
	PHDFilter
include("classes.jl")

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

end # module
