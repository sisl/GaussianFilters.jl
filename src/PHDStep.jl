__precompile__(true)
module PHDFiltering

export step, prune, step_prune



"""
Function to combine stepping through PHD filter and pruning resultant distribution

Arguments:
- `x::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Matrix` Matrix of measurements. [Measurements]
- `PHD::PHDFilter` PHD filter to step through. [PHD Filter]
- `T::Real` Threshholding value for pruning. [probability]

Returns:
- `xp::GaussianMixture` Describe first return value. [Gaussian Mixture]

Notes:
1. Make specific comments about corner cases or things to know about usage here.

References:
1. List citations to references used in writing the function here
"""
function step(x::GaussianMixture, Z::Matrix, PHD::PHDFilter)
    
    J_k = x.N
    _, m = size(x.w)
    w_p = zeros(0)
    μ_p = zeros(0,m)
    Σ_p = zeros(0,m,m)
    
    # add expand along first dimension
    add_dim0(x::Array) = reshape(x, (1,size(x)...))

    # 1. Prediction for birth targets
    J_γ = PHD.γ.N
    for j in 1:J_γ
       w_p = [w_p; PHD.γ.w[j]]
       μ_p = cat(μ_p, add_dim0(PHD.γ.μ[j,:]); dims=1)
       Σ_p = cat(Σ_p, add_dim0(PHD.γ.Σ[j,:,:]); dims=1)   
    end

    J_β = PHD.spawn.β.N
    for j in 1:J_β
        for l in 1:J_k
            spawn_dyn = PHD.spawn.dyn[j]
            w_p = [w_p x.w[l]*PHD.spawn.β.w[j]]

            μ_add = spawn_dyn.d + spawn_dyn.A*x.μ[l,:]
            μ_p = cat(μ_p, add_dim0(μ_add); dims=1)

            Σ_add = spawn_dyn.Q + spawn_dyn.A*x.Σ[l,:,:]*spawn_dyn.A'
            Σ_p = cat(Σ_p, add_dim0(Σ_add); dims=1) 
        end
    end

    # 2. Prediction for existing targets
    for j in 1:J_k
        w_p = [w_p PHD.Ps*x.w[j]]

        μ_p = cat(μ_p, add_dim0(PHD.dyn.A*x.μ[j,:]); dims=1)

        Σ_add = PHD.dyn.Q + PHD.dyn.A*x.Σ[j,:,:]*PHD.dyn.A'
        Σ_p = cat(Σ_p, add_dim0(Σ_add); dims=1) 
    end

    x_p = GaussianMixture(w_p, μ_p, Σ_p)
    # 3. Construction of PHD update components
    

    # 4. Update
      
    
    
    x_n = GaussianMixture(w_n, μ_n, Σ_n)
    return x_n
end




"""
INSERT PRUNE
"""


"""
Function to combine stepping through PHD filter and pruning resultant distribution

Arguments:
- `x::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Matrix` Matrix of measurements. [Measurements]
- `PHD::PHDFilter` PHD filter to step through. [PHD Filter]
- `T::Real` Threshholding value for pruning. [probability]

Returns:
- `xp::GaussianMixture` Describe first return value. [Gaussian Mixture]

Notes:
1. Make specific comments about corner cases or things to know about usage here.

References:
1. List citations to references used in writing the function here
"""
function step_prune(x::GaussianMixture, Z::Matrix, PHD::PHDFilter, T::Real)
    xp = step(x, Z, PHD)
    xpp = prune(xp, T)
    return xp
end






end # End of module PHDStep