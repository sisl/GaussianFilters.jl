# Implementation of the Interacting Multiple Model filter 
# Reference: Bar-Shalom, Yaakov, X. Rong Li, and Thiagalingam Kirubarajan. 
# Estimation with applications to tracking and navigation: theory algorithms and software. 
# John Wiley & Sons, 2004.

struct InteractingMultipleModelFilter{F<:AbstractFilter, D<:DynamicsModel, O<:ObservationModel, T<:AbstractMatrix} <: AbstractFilter
    filters::Vector{F}
    dynamics_models::Vector{D}
    observation_models::Vector{O}
    transition::T
end

struct InteractingMultipleModelBelief{GB<:GaussianBelief, T<:Real}
    beliefs::Vector{GB}
    models_probs::Vector{T}
end

function update(filter::InteractingMultipleModelFilter, b0::InteractingMultipleModelBelief, 
                u::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    # predict
    bp = predict(filter, b0, u)

    # measure
    bn = measure(filter, bp, y; u = u)
end

function predict(filter::InteractingMultipleModelFilter, b0::InteractingMultipleModelBelief{GB,T},
                 u::AbstractVector{<:Real}) where {GB<:GaussianBelief, T<:Real}
    
    nm = length(b0.beliefs)
    new_beliefs = Vector{GB}(undef, nm)
    for i=1:nm
        new_beliefs[i] = predict(filter.filters[i], b0.beliefs[i], u)
    end

end

function mixing(filter::InteractingMultipleModelFilter, b0::InteractingMultipleModelBelief{GB,T},
                 u::AbstractVector{<:Real}) where {GB<:GaussianBelief, T<:Real}
                 
    # mixing probabilities
    mixing_probs = filter.transition .* b0.models_probs
    mixing_probs = mixing_probs ./ sum(mixing_probs, dims=1)

    # mixing
    nm = length(b0.beliefs)
    mixed_beliefs = Vector{GB}(undef, nm)
    for j=1:nm
        mixed_mu = sum(mixing_probs[i,j] .* b0.beliefs[i].μ  for i=1:nm)
        mixed_cov = sum(
            mixing_probs[i,j].*(b0.beliefs[i].Σ + 
            (b0.beliefs[i].μ - mixed_mu[j])*(b0.beliefs[i].μ - mixed_mu[j])')
            for i=1:nm
        )
        mixed_beliefs[j] = GaussianBelief(mixed_mu, mixed_cov)
    end
    return mixed_beliefs
end

function mode_filtering(filter::InteractingMultipleModelFilter,
                 b0::InteractingMultipleModelBelief{GB,T},
                 u::AbstractVector{<:Number}, 
                 y::AbstractVector{<:Number}) where {GB<:GaussianBelief, T<:Number}
    nm = length(b0.beliefs)
    new_beliefs = Vector{GB}(undef, nm)
    for j=1:nm
        new_beliefs[j] = update(filter, b0, u, y)
    end
end

function mode_likelihood(filter::InteractingMultipleModelFilter, b0::InteractingMultipleModelBelief{GB,T},
                 u::AbstractVector{<:Number}, y::AbstractVector{<:Number}) where {GB<:GaussianBelief, T<:Number}
    nm = length(b0.beliefs)
    for j=1:nm
    end
end

function mode_probability_update(filter::InteractingMultipleModelFilter,
                 b0::InteractingMultipleModelBelief{GB,T},
                 u::AbstractVector{<:Number}, 
                 y::AbstractVector{<:Number}) where {GB<:GaussianBelief, T<:Number}
    nm = length(b0.beliefs)
    for j=1:nm
        # measurement prediction 
        z = measure(filter.models[j], )
    end
end