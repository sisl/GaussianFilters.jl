__precompile__(true)
module GaussianFilters

# Usings
using Reexport,
LinearAlgebra

# Base Module Includes
include("classes.jl")
include("PHDStep.jl")

# Export Module Contents
#@reexport using GaussianFilters.YourSubmodule
#@reexport using GaussianFilters.RubberDucks
@reexport using GaussianFilters.Classes
@reexport using GaussianFilters.PHDFilter


end # module
