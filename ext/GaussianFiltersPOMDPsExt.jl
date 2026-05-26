module GaussianFiltersPOMDPsExt

using GaussianFilters
using POMDPs

# A GaussianFilters AbstractFilter doubles as a POMDPs.Updater. The two
# update signatures already agree on argument order: (b_prev, a, o).
POMDPs.update(f::AbstractFilter, b::GaussianBelief, a::AbstractVector, o::AbstractVector) =
    GaussianFilters.update(f, b, a, o)

# Pass-through initialize_belief for callers who hold a GaussianBelief.
POMDPs.initialize_belief(::AbstractFilter, b::GaussianBelief) = b

end # module
