### Measurement Model ###
"""
    Measurement(C,R)

Construct measurement model with observation matrix C and sensor
noise matrix R
"""
mutable struct Measurement{a,b}
    C::Matrix{a}
    R::Matrix{b}
end

### Dynamics Model ###
"""
    Dynamics(A,Q,d)
    Dynamics(A,Q)

Construct linear dynamics model with; transition matrix A,
process noise matrix Q and offset vector d
"""
mutable struct Dynamics{a,b,c}
    A::Matrix{a}
    Q::Matrix{b}
    d::Vector{c}
end

## Constructors ##
function Dynamics(A,Q)
    return Dynamics(A,Q,Int8[])
end

"""
    GaussianMixture(N, w, μ, Σ)

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
    @assert length(μ) == length(Σ) == length(w) "bad length"
    GaussianMixture{promote_type(T,K)}(length(w), w, μ, Σ)
end

"""
    Spawn(β, dyn)

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
    PHDFilter(γ, spawn, dyn, meas, Ps, Pd)

    Sets up a PHD Filter

    Arguments:
    γ: Birth intensity
    spawn: Spawning intensity
    dyn: Dynamics
    meas: Measurements
    Ps: Survival probability
    Pd: Detection probability
"""

struct PHDFilter
    γ::GaussianMixture
    spawn::Spawn
    dyn::Dynamics
    meas::Measurement
    Ps::Float64
    Pd::Float64
end
