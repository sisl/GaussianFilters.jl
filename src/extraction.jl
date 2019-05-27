"""
    Multiple-Target State Extraction(x, threshold)

    Extracts targets whose weights (x.w) are above threshold.

    Arguments:
        x: Set of Gaussian Mixtures
        threshold: Threshold on weights. Above this threshold,
        state estimate is extracted

    Returns:
        x_new: Extracted set of Gaussian Mixtures
"""

function multiple_target_state_extraction(x, threshold)
    
    inds = x.w .> threshold
    
    return GaussianMixture(x.w[inds], x.μ[inds], x.Σ[inds])        
end
