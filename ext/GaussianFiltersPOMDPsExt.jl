module GaussianFiltersPOMDPsExt

using GaussianFilters
using POMDPs
using Distributions: AbstractMvNormal, mean, cov

# Two integration patterns are supported:
#
# 1) Direct dispatch on AbstractFilter — convenient when you just want
#    POMDPs.update(filter, b, a, o) without going through the full
#    POMDPs.simulate machinery.
#
# 2) A GaussianFilterUpdater wrapper subtyping POMDPs.Updater — needed
#    by simulators (HistoryRecorder, RolloutSimulator) that dispatch on
#    the Updater abstract type.

# --- direct AbstractFilter methods ---

POMDPs.update(f::AbstractFilter, b::GaussianBelief, a::AbstractVector, o::AbstractVector) =
    GaussianFilters.update(f, b, a, o)

POMDPs.initialize_belief(::AbstractFilter, b::GaussianBelief) = b
POMDPs.initialize_belief(::AbstractFilter, d::AbstractMvNormal) =
    GaussianBelief(mean(d), cov(d))

# --- POMDPs.Updater wrapper ---

"""
    GaussianFilterUpdater(filter::AbstractFilter)

A `POMDPs.Updater` wrapper around any GaussianFilters `AbstractFilter`.
Use this when handing a Gaussian-filter belief updater to a POMDPs.jl
simulator that requires `<:Updater` dispatch
(`HistoryRecorder`, `RolloutSimulator`, etc).
"""
struct GaussianFilterUpdater{F<:AbstractFilter} <: POMDPs.Updater
    filter::F
end

POMDPs.update(u::GaussianFilterUpdater, b::GaussianBelief, a::AbstractVector, o::AbstractVector) =
    GaussianFilters.update(u.filter, b, a, o)

POMDPs.initialize_belief(::GaussianFilterUpdater, b::GaussianBelief) = b
POMDPs.initialize_belief(::GaussianFilterUpdater, d::AbstractMvNormal) =
    GaussianBelief(mean(d), cov(d))

# Make the wrapper reachable from a `using GaussianFilters` context by
# overloading a stub function defined in the main package.
GaussianFilters.pomdps_updater(f::AbstractFilter) = GaussianFilterUpdater(f)

end # module
