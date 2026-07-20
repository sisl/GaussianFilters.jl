### Filters ###

"""
`AbstractFilter` is an abstract type to encapsulate different kinds
of discrete gaussian filters
"""
abstract type AbstractFilter end


"""
    KalmanFilter(d::LinearDynamicsModel,o::LinearObservationModel)

Construct Kalman filter with LinearDynamicsModel d and
LinearObservationModel o.
"""
struct KalmanFilter{D<:LinearDynamicsModel, O<:LinearObservationModel} <: AbstractFilter
    d::D
    o::O
end

"""
    ExtendedKalmanFilter(d::DynamicsModel,o::ObservationModel)
    KalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number,
        α::Float,β::Float)

Construct Extended Kalman filter with DynamicsModel d and
ObservationModel o.
"""
struct ExtendedKalmanFilter{D<:DynamicsModel, O<:ObservationModel} <: AbstractFilter
    d::D
    o::O
    function ExtendedKalmanFilter{D,O}(d::D, o::O) where {D<:DynamicsModel, O<:ObservationModel}
        (d isa LinearDynamicsModel && o isa LinearObservationModel) &&
            error("linear models: use kf over ekf")
        return new{D,O}(d, o)
    end
end

ExtendedKalmanFilter(d::D, o::O) where {D<:DynamicsModel, O<:ObservationModel} =
    ExtendedKalmanFilter{D,O}(d, o)

"""
    UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number,
        α::Float,β::Float)
    UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number)
    UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel)

Construct Unscented Kalman filter with DynamicsModel d, ObservationModel o,
and UKF parameters λ, α, and β. Default constructor uses α/β formulation from
Probabilistic Robotics, second constructor reduces complexity, third
constructor defaults λ to 2, as is commonly done.
"""
struct UnscentedKalmanFilter{D<:DynamicsModel, O<:ObservationModel,
                              a<:Number, b<:Number, c<:Number} <: AbstractFilter
    d::D
    o::O
    λ::a
    α::b
    β::c
end

function UnscentedKalmanFilter(d::DynamicsModel, o::ObservationModel, λ::Number)
    return UnscentedKalmanFilter(d, o, λ, 1, 0)
end

function UnscentedKalmanFilter(d::DynamicsModel, o::ObservationModel)
    return UnscentedKalmanFilter(d, o, 2, 1, 0)
end

### Belief States ###

"""
    GaussianBelief(μ::AbstractVector,Σ::Symmetric)
    GaussianBelief(μ::AbstractVector,Σ::AbstractMatrix)

Construct a gaussian belief, consisting of mean vector μ
and symmetric covariance matrix Σ
"""
struct GaussianBelief{T<:AbstractVector{<:Number}, S<:Symmetric{<:Number}}
    μ::T
    Σ::S
    #=
    function GaussianBelief{a,b}(μ::Vector{a},
        Σ::Symmetric{b}) where {a<:Number,b<:Number}
        (size(μ,1) == size(Σ,1)) ?
        new(μ,Σ) : error("first dimensions not matching")
    end
    =#
end

function GaussianBelief(μ::AbstractVector,Σ::AbstractMatrix)
    return GaussianBelief(μ,Symmetric(Σ))
end

function GaussianBelief(μ::AbstractVector, Σ::UniformScaling)
    return GaussianBelief(μ, Σ(length(μ)))
end

function Base.rand(rng::AbstractRNG, b::GaussianBelief)
    return b.μ + cholesky(b.Σ).L * randn(rng, size(b.Σ,1))
end
