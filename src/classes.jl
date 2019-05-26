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
