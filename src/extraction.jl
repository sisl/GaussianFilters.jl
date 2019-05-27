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

function MultipleTargetStateExtraction(x, threshold)::GaussianMixture

    Jk = length(x.w)    
    new_w = zeros(Jk)
    new_μ = zeros(size(x.μ))
    new_Σ = zeros(size(x.Σ))
    
    j = 1
    for i=1:Jk
        if x.w[i] > threshold
            new_w[j] = w[i]
            new_μ[j] = μ[i]
            new_Σ[j,:,:] = Σ[i,:,:]
            j += 1
        end
    end
    
    return GaussianMixture(new_w,new_μ,new_Σ)        
end
