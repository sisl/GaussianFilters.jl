"""
    GaussianMixture(N, w, μ, Σ)

    Arguments:
    N: number of models
    μ: means of the model
    Σ: Covariances of the model
"""

struct GaussianMixture
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

