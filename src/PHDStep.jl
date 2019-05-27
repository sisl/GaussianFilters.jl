__precompile__(true)
module PHDFiltering

export step, prune, step_prune

"""
Function to evaluate multivariate normal pdf function parameterized by μ
and Σ at x

Arguments:
- `x::Vector` Point to evaluate pdf
- `μ::Vector` Mean of Multivariate Normal Gaussian
- `Σ::Matrix` Covariance Multivariate Normal Gaussian

Returns:
- `xp::Vector` evaluation of multivariate normal pdf
"""
MvNormalPDF(x::Vector, μ::Vector, Σ::Matrix) = (2*pi)^(-length(μ)/2)*det(Σ)^(-1/2)*exp(-1/2*(x-μ)'*inv(Σ)*(x-μ))


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
function step(x::GaussianMixture, Z::Vector{Vector}, PHD::PHDFilter)

    J_k = x.N
    w_p = []
    μ_p = []
    Σ_p = []
    
    # 1. Prediction for birth targets
    J_γ = PHD.γ.N
    for j in 1:J_γ
        push!(w_p,PHD.γ.w[j])
        push!(μ_p,PHD.γ.μ[j])
        push!(Σ_p,PHD.γ.Σ[j]) 
    end

    J_β = PHD.spawn.β.N
    for j in 1:J_β, ℓ in 1:J_k
        spawn_dyn = PHD.spawn.dyn[j]
        push!(w_p, x.w[ℓ]*PHD.spawn.β.w[j])

        μ_add = spawn_dyn.d + spawn_dyn.A*x.μ[ℓ]
        push(μ_p, μ_add)

        Σ_add = spawn_dyn.Q + spawn_dyn.A*x.Σ[ℓ]*spawn_dyn.A'
        push!(Σ_p, Σ_add) 
    end

    # 2. Prediction for existing targets
    for j in 1:J_k
        push!(w_p, PHD.Ps*x.w[j])

        push!(μ_p, PHD.dyn.A*x.μ[j])

        Σ_add = PHD.dyn.Q + PHD.dyn.A*x.Σ[j]*PHD.dyn.A'
        push!(Σ_p, Σ_add) 
    end
    x_p = GaussianMixture(w_p, μ_p, Σ_p)

    # 3. Construction of PHD update components
    
    J_p = x_p.N

    η_n = []
    S_n = []
    K_n = []
    P_n = []
    for j in 1:J_p
        push!(η_n, PHD.meas.C*x_p.μ[j])

        S_add = PHD.meas.R + PHD.meas.C*x_p.Σ[j]*PHD.meas.C'
        push!(S_n, S_add) 

        K_add = x_p.Σ[j]*PHD.meas.C'*inv(S_add)
        push!(K_n, K_add)

        P_add = (1.0I - K_add*PHD.meas.C)*x_p.Σ[j]
        push!(P_n, P_add) 
    end

    # 4. Update
    w_n = []
    μ_n = []
    Σ_n = []
    for j in 1:J_p
        push!(w_n,(1-Pd)*w_p[j])
        push!(μ_n, x_p.μ[j])
        push!(Σ_n, x_p.Σ[j])
    end

    for z in Z
        w_ntemp = []
        for j in 1:J_p
            push!(w_ntemp, Pd*x_p.w[j]*MvNormalPDF(z,η_n[j],Σ_n[j]))
            push!(μ_n, x_p[j] + K_n[j]*(z-η_n[j]))
            push!(Σ_n, P_n[j])
        end
        push!(w_n, w_ntemp ./ (PDF.κ(z) + sum(w_ntemp)))
    end  
  
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