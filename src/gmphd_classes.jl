"""
    Measurement(C::AbstractMatrix, R::AbstractMatrix)

Construct measurement model with observation matrix C and sensor
noise matrix R
"""
struct Measurement{A<:Number,B<:Number}
    C::AbstractMatrix{A}
    R::AbstractMatrix{B}
end

"""
    Dynamics(A::AbstractMatrix, Q::AbstractMatrix, d::AbstractVector)
    Dynamics(A::AbstractMatrix, Q::AbstractMatrix)

Construct linear dynamics model with; transition matrix A,
process noise matrix Q and offset vector d
"""
struct Dynamics{A<:Number,B<:Number,C<:Number}
    A::AbstractMatrix{A}
    Q::AbstractMatrix{B}
    d::AbstractVector{C}
end

function Dynamics(A,Q)
    N = size(Q,1)
    return Dynamics(A,Q,zeros(N))
end

"""
    GaussianMixture(N::Int64, w::Vector{Number}, μ::Vector{AbstractVector},
        Σ::Vector{AbstractMatrix})
    GaussianMixture(w::Vector{Number}, μ::Vector{AbstractVector},
        Σ::Vector{AbstractMatrix})

    Arguments:
    N: Number of models
    w: Weights
    μ: Means of the model
    Σ: Covariances of the model
"""
struct GaussianMixture{T<:Number}
    N::Int64
    w::Vector{<:Number}
    μ::Vector{<:AbstractVector{T}}
    Σ::Vector{<:AbstractMatrix{T}}
end

function GaussianMixture(w, μ::Vector{<:AbstractVector{T}},
    Σ::Vector{<:AbstractMatrix{K}}) where {T<:Number,K<:Number}
    @assert length(μ) == length(Σ) == length(w) "Number
     of mixtures inconsistent"
    GaussianMixture{promote_type(T,K)}(length(w), w, μ, Σ)
end

"""
    Spawn(β::GaussianMixture, dyn::Vector{Dynamics})

    Arguments:
    β: Gaussian Mixture model determining the
    spawning intesity of the target
    dyn: Dynamics
"""
struct Spawn
    β::GaussianMixture
    dyn::Vector{Dynamics}
end

"""
    PHDFilter(γ::GaussianMixture, spawn::Spawn, dyn::Vector{Dynamics},
        meas::Measurement, Ps::Float64, Pd::Float64, κ::Function)
    PHDFilter(γ::GaussianMixture, spawn::Spawn, dyn::Dynamics,
        meas::Measurement, Ps::Float64, Pd::Float64, κ::Function)

    Sets up a PHD Filter

    Arguments:
    γ: Birth intensity
    spawn: Spawning intensity
    dyn: Vector of possible Dynamics models
    meas: Measurements
    Ps: Survival probability
    Pd: Detection probability
"""
struct PHDFilter
    γ::GaussianMixture
    spawn::Spawn
    dyn::Vector{Dynamics}
    meas::Measurement
    Ps::Float64
    Pd::Float64
    κ::Function
end

function PHDFilter(γ::GaussianMixture, spawn::Spawn, dyn::Dynamics,
    meas::Measurement, Ps::Float64, Pd::Float64, κ::Function)
    return PHDFilter(γ, spawn, [dyn], meas, Ps, Pd, κ)
end
