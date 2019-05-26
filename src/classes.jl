### Dynamics Model ###
"""
    Dynamics(A,Q,d)
    Dynamics(A,Q)

Construct linear dynamics model with; transition matrix A,
process noise matrix Q and constant matrix d
"""
mutable struct Dynamics
    A::Matrix{Float64}
    Q::Matrix{Float64}
    d::Vector{Float64}
end

## Constructors ##
function Dynamics(A,Q)
    return Dynamics(A,Q,d = nothing)
end

function Dynamics(A,Q,d)
    return Dynamics(A,Q,d)
end

### Measurement Model ###
"""
    Measurement(C,R)

Construct measurement model with observation matrix C and sensor
noise matrix R
"""
mutable struct Measurement
    C::Matrix{Float64}
    R::Matrix{Float64}
end

## Constructor ##
function Measurement(C,R)
    return Measurement(C,R)
end

"""
    GaussianMixture(N, w, μ, Σ)

    Arguments:
    N: Number of models
    w: Weights
    μ: Means of the model
    Σ: Covariances of the model
"""

struct GaussianMixture
    N :: Int64
    w :: Vector{Float64}
    μ :: Matrix{Float64}
    Σ :: Array{Float64,3}
end

function GaussianMixture(w,μ,Σ)
    return GaussianMixture(N = length(w), w, μ, Σ)
end

function GaussianMixture(N,w,μ,Σ)
    return GaussianMixture(N,w,μ,Σ)
end

"""
    Spawn(β, dyn)

    Arguments:
    β: Gaussian Mixture model determining the 
    spawning intesity of the target
    dyn: Dynamics
"""

struct Spawn
    β :: GaussianMixture
    dyn :: Dynamics
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
    γ :: GaussianMixture
    spawn :: Spawn
    dyn :: Dynamics
    meas :: Measurement
    Ps :: Float64
    Pd :: Float64
end
