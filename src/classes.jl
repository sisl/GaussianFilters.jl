"""
    GaussianMixture(N, w, μ, Σ)

    Arguments:
    N: number of models
    μ: means of the model
    Σ: Covariances of the model
"""

struct GaussianMixture{N, w, μ, Σ}
    N :: Int64
    w :: Vector{Float64}
    μ :: Matrix{Float64,2}
    Σ :: Array{Float64,3}
end

"""
    Spawn(β, dyn)

    Arguments:
    β: Gaussian Mixture model determining the
    spawning intesity of the target
    dyn: Dynamics
"""

struct Spawn{β, dyn}
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

struct PHDFilter{γ, spawn, dyn, meas, Ps, Pd}
    γ :: GaussianMixture
    spawn :: Spawn
    dyn :: Dynamics
    meas :: Measurement
    Ps :: Float64
    Pd :: Float64
end

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
