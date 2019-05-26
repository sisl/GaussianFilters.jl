"""
    GaussianMixture

"""

struct GaussianMixture{N,ω,μ,Σ}
    N :: Int64
    ω :: Vector{Float64}
    μ :: Matrix{Float64,2}
    Σ :: Array{Float64,3}
end
