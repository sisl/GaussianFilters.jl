"""
    Multiple-Target State Extraction(x, threshold)

    Extracts targets whose weights (x.w) are above threshold.

    Arguments:
        x: Set of Gaussian Mixtures
        threshold: Threshold on weights. Above this threshold,
        state estimate is extracted

    Returns:
        X: Multi-Target State Estimate
"""

function multiple_target_state_extraction(x, threshold)
    inds = x.w .> threshold
    N = sum(inds)
    m = x.μ[inds]
    w = x.w[inds]

    Xk = []
    for i = 1:N
        tmp = [m[i] for j = 1:round(w[i])]
        Xk = vcat(Xk, tmp)
    end

    return Xk
end
