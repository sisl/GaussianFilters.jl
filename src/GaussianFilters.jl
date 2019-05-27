__precompile__(true)
module GaussianFilters

# Usings
using Reexport
using LinearAlgebra
using Random
import Random: rand

export
	Measurement,
	Dynamics,
	GaussianMixture,
	Spawn,
	PHDFilter,
	step,
	prune,
	step_prune,
	multiple_target_state_extraction,

# Base Module Includes
include("classes.jl")
include("PHDStep.jl")
include("prune.jl")
include("extraction.jl")

# Export Module Contents
#@reexport using GaussianFilters.YourSubmodule
#@reexport using GaussianFilters.RubberDucks
#@reexport using GaussianFilters.Classes
#@reexport using GaussianFilters.PHDFilter


end # module
