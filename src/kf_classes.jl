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
struct KalmanFilter <: AbstractFilter
    d::LinearDynamicsModel
    o::LinearObservationModel
end

"""
    ExtendedKalmanFilter(d::DynamicsModel,o::ObservationModel)
    KalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number,
        α::Float,β::Float)

Construct Extended Kalman filter with DynamicsModel d and
ObservationModel o.
"""
struct ExtendedKalmanFilter <: AbstractFilter
    d::DynamicsModel
    o::ObservationModel
    ExtendedKalmanFilter(d,o) =
        (d isa LinearDynamicsModel && o isa LinearObservationModel) ?
        error("linear models: use kf over ekf") : new(d,o)
end

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
struct UnscentedKalmanFilter{a<:Number,b<:Number,c<:Number} <: AbstractFilter
    d::DynamicsModel
    o::ObservationModel
    λ::a
    α::b
    β::c
    #=
    UnscentedKalmanFilter(d,o,λ,α,β) = (d isa LinearDynamicsModel && o isa LinearObservationModel) ?
        error("linear models: use kf over ukf") : new(d,o,λ,α,β)
    =#
end

function UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::a) where a<:Number
    return UnscentedKalmanFilter{a,Int8,Int8}(d,o,λ,1,0)
end

function UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel)
    return UnscentedKalmanFilter{Int8,Int8,Int8}(d,o,2,1,0)
end

### Belief States ###

"""
    GaussianBelief(μ::AbstractVector,Σ::Symmetric)
    GaussianBelief(μ::AbstractVector,Σ::AbstractMatrix)

Construct a gaussian belief, consisting of mean vector μ
and symmetric covariance matrix Σ
"""
struct GaussianBelief{a<:Number,b<:Number}
    μ::AbstractVector{a}
    Σ::Symmetric{b}
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

function Base.rand(rng::AbstractRNG, b::GaussianBelief)
    return b.μ + cholesky(b.Σ).L * randn(size(b.Σ,1))
end
