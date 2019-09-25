"""
    multiple_target_state_extraction(b::GaussianMixture, th::Real)

Extracts targets whose weights (x.w) are above threshold.

Arguments:
    `b::GaussianMixture`: Set of Gaussian Mixtures
    `th::Real`: Threshold on weights. Above this threshold,
    state estimate is extracted

Returns:
    X: Multi-Target State Estimate
"""
function multiple_target_state_extraction(b::GaussianMixture, th::Real)
    inds = b.w .> th
    N = sum(inds)
    m = b.μ[inds]
    w = b.w[inds]

    X = []
    for i = 1:N
        tmp = [m[i] for j = 1:round(w[i])]
        X = vcat(X, tmp)
    end

    return X
end
