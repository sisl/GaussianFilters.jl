"""
    GaussianMixture(N,w,μ,Σ)

    Arguments:
    N: number of models
    μ: means of the model
    Σ: Covariances of the model
"""

struct GaussianMixture{N,w,μ,Σ}
    N :: Int64
    w :: Vector{Float64}
    μ :: Matrix{Float64,2}
    Σ :: Array{Float64,3}
end

