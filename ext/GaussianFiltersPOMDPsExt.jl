module GaussianFiltersPOMDPsExt

using GaussianFilters
using POMDPs
using Distributions: AbstractMvNormal, mean, cov

# A GaussianFilters AbstractFilter doubles as a POMDPs.Updater. The two
# update signatures already agree on argument order: (b_prev, a, o).
POMDPs.update(f::AbstractFilter, b::GaussianBelief, a::AbstractVector, o::AbstractVector) =
    GaussianFilters.update(f, b, a, o)

# Initialize from a GaussianBelief (identity) or from any multivariate
# normal distribution (extract mean and covariance).
POMDPs.initialize_belief(::AbstractFilter, b::GaussianBelief) = b
POMDPs.initialize_belief(::AbstractFilter, d::AbstractMvNormal) =
    GaussianBelief(mean(d), cov(d))

end # module
