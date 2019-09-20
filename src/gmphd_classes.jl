"""
    Measurement(C::Matrix, R::Matrix)

Construct measurement model with observation matrix C and sensor
noise matrix R
"""
mutable struct Measurement{a,b}
    C::Matrix{a}
    R::Matrix{b}
end

"""
    Dynamics(A::Matrix, Q::Matrix, d::Vector)
    Dynamics(A::Matrix, Q::Matrix)

Construct linear dynamics model with; transition matrix A,
process noise matrix Q and offset vector d
"""
mutable struct Dynamics{a,b,c}
    A::Matrix{a}
    Q::Matrix{b}
    d::Vector{c}
end

function Dynamics(A,Q)
    N = size(Q,1)
    return Dynamics(A,Q,zeros(N))
end

"""
    GaussianMixture(N::Int64, w::Vector{Float64}, μ::Vector{Vector},
        Σ::Vector{Matrix})
    GaussianMixture(w::Vector{Float64}, μ::Vector{Vector}, Σ::Vector{Matrix})

    Arguments:
    N: Number of models
    w: Weights
    μ: Means of the model
    Σ: Covariances of the model
"""
struct GaussianMixture{T<:Number}
    N::Int64
    w::Vector{Float64}
    μ::Vector{Vector{T}}
    Σ::Vector{Matrix{T}}
end

function GaussianMixture(w, μ::Vector{Vector{T}}, Σ::Vector{Matrix{K}}) where {T,K}
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
