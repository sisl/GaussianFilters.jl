__precompile__(true)
module JuliaPackageTemplate

# Usings
using Reexport

# Base Module Includes
include("submodule.jl")
include("rubber_ducks.jl")

# Export Module Contents
@reexport using JuliaPackageTemplate.YourSubmodule
@reexport using JuliaPackageTemplate.RubberDucks

end # module
