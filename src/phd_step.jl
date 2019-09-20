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

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function step(x::GaussianMixture, Z::Vector{Vector{T}},
                PHD::PHDFilter) where T <: Real

    J_k = x.N
    n = length(x.μ[1])
    J_γ = PHD.γ.N
    J_β = PHD.spawn.β.N
    J_dyn = length(PHD.dyn)
    tot = J_γ + J_β*J_k + J_k*J_dyn
    w_p = [0.0 for i in 1:tot]
    μ_p = [zeros(n) for i in 1:tot]
    Σ_p = [zeros(n,n) for i in 1:tot]

    # 1. Prediction for birth targets
    i = 0
    for j in 1:J_γ
        i += 1
        w_p[i] = PHD.γ.w[j]
        μ_p[i] = PHD.γ.μ[j]
        Σ_p[i] = PHD.γ.Σ[j]
    end


    for j in 1:J_β, ℓ in 1:J_k
        i += 1
        spawn_dyn = PHD.spawn.dyn[j]
        w_p[i] = x.w[ℓ]*PHD.spawn.β.w[j]
        μ_p[i] =  spawn_dyn.d + spawn_dyn.A*x.μ[ℓ]
        Σ_p[i] = spawn_dyn.Q + spawn_dyn.A*x.Σ[ℓ]*spawn_dyn.A'
    end

    # 2. Prediction for existing targets
    for j in 1:J_k
        for k in 1:J_dyn
            i += 1
            w_p[i] = PHD.Ps*x.w[j]
            μ_p[i] = PHD.dyn[k].A*x.μ[j]
            Σ_p[i] = PHD.dyn[k].Q + PHD.dyn[k].A*x.Σ[j]*PHD.dyn[k].A'
        end
    end
    @assert i == tot "Poor Indexing"
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
    w_n = similar(x_p.w,0)
    μ_n = similar(x_p.μ,0)
    Σ_n = similar(x_p.Σ,0)
    for j in 1:J_p
        push!(w_n,(1-PHD.Pd)*w_p[j])
        push!(μ_n, x_p.μ[j])
        push!(Σ_n, x_p.Σ[j])
    end

    for z in Z
        w_ntemp = []
        for j in 1:J_p
            push!(w_ntemp, PHD.Pd*x_p.w[j]*MvNormalPDF(z,η_n[j],S_n[j]))
            push!(μ_n, x_p.μ[j] + K_n[j]*(z-η_n[j]))
            push!(Σ_n, P_n[j])
        end
        append!(w_n, w_ntemp ./ (PHD.κ(z) + sum(w_ntemp)))
    end

    x_n = GaussianMixture(w_n, μ_n, Σ_n)
    return x_n
end

"""
Function to combine stepping through PHD filter and pruning resultant distribution

Arguments:
- `x::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Matrix` Matrix of measurements. [Measurements]
- `PHD::PHDFilter` PHD filter to step through. [PHD Filter]
- `T::Real` Threshholding value for pruning. [probability]

Returns:
- `xp::GaussianMixture` Describe first return value. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function step_prune(x::GaussianMixture, Z::Vector{Vector{N}}, PHD::PHDFilter,
    T::Real, U::Real, J_max::Integer) where N <: Real
    xp = step(x, Z, PHD)
    xpp = prune(xp, T, U, J_max)
    return xpp
end
