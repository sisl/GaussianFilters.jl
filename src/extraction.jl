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

function MultipleTargetStateExtraction(x, threshold)
    
    Jk = length(x.w)    
    new_w = []
    new_μ = []
    new_Σ = []

    for i=1:Jk
        if x.w[i] > threshold
            push!(new_w, w[i])
            push!(new_μ, μ[i])
            push!(new_Σ, Σ[i])
        end
    end
    
    return GaussianMixture(new_w,new_μ,new_Σ)        
end
