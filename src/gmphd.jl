"""
	update(x::GaussianMixture, Z::Vector{AbstractVector}, PHD::PHDFilter,
		T::Real, U::Real, J_max::Integer)

Perform an update step using a GMPHD Filter.

Arguments:
- `x::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Matrix` Matrix of measurements. [Measurements]
- `PHD::PHDFilter` PHD filter to step through. [PHD Filter]
- `T::Real` Truncation threshold.  Drop distributions with weight less than T
- `U::Real` Merging threshold.  Merge if (μ_1-μ_2)^T * Σ^-1 * (μ_1-μ_2) < U
- `J_max::Integer` Maximum number of features

Returns:
- `xn::GaussianMixture` Describe first return value. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function update(x::GaussianMixture, Z::Vector{<:AbstractVector{<:Real}},
	PHD::PHDFilter, T::Real, U::Real, J_max::Integer)

	xp = predict(x, PHD)
	xm = measure(xp, Z, PHD)
    xn = prune(xm, T, U, J_max)
    return xn
end

"""
	predict(x::GaussianMixture, PHD::PHDFilter)

Make a prediction on next state based on PHD dynamics

Arguments:
- `x::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `PHD::PHDFilter` PHD filter to step through. [PHD Filter]

Returns:
- `xp::GaussianMixture` Predicted next state distribution. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function predict(x::GaussianMixture, PHD::PHDFilter)

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
    xp = GaussianMixture(w_p, μ_p, Σ_p)
	return xp
end

"""
	measure(xp::GaussianMixture, Z::Vector{AbstractVector}, PHD::PHDFilter)

Perform a measurement update on a predicted PHD next state.

Arguments:
- `xp::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Vector{AbstractVector}` Array of measurements. [Measurements]
- `PHD::PHDFilter` PHD filter to step through. [PHD Filter]

Returns:
- `xm::GaussianMixture` Describe first return value. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function measure(xp::GaussianMixture, Z::Vector{<:AbstractVector{T}},
                PHD::PHDFilter) where T <: Real

	# 3. Construction of PHD update components

    J_p = xp.N

    η_n = []
    S_n = []
    K_n = []
    P_n = []
    for j in 1:J_p
        push!(η_n, PHD.meas.C*xp.μ[j])

        S_add = PHD.meas.R + PHD.meas.C*xp.Σ[j]*PHD.meas.C'
        push!(S_n, S_add)

        K_add = xp.Σ[j]*PHD.meas.C'*inv(S_add)
        push!(K_n, K_add)

        P_add = (1.0I - K_add*PHD.meas.C)*xp.Σ[j]
        push!(P_n, P_add)
    end

    # 4. Update
    w_n = similar(xp.w,0)
    μ_n = similar(xp.μ,0)
    Σ_n = similar(xp.Σ,0)
    for j in 1:J_p
        push!(w_n,(1-PHD.Pd)*xp.w[j])
        push!(μ_n, xp.μ[j])
        push!(Σ_n, xp.Σ[j])
    end

    for z in Z
        w_ntemp = []
        for j in 1:J_p
            push!(w_ntemp, PHD.Pd*xp.w[j]*MvNormalPDF(z,η_n[j],S_n[j]))
            push!(μ_n, xp.μ[j] + K_n[j]*(z-η_n[j]))
            push!(Σ_n, P_n[j])
        end
        append!(w_n, w_ntemp ./ (PHD.κ(z) + sum(w_ntemp)))
    end

    xm = GaussianMixture(w_n, μ_n, Σ_n)
    return xm
end

"""
    prune(x::GaussianMixture, T::Real, U::Real, J_max::Integer)

Prune the posterior Gaussian-Mixture next PHD state

Arguments:
- `x::GaussianMixture` Posterior state distribution. [Gaussian Mixture]
- `T::Real` Truncation threshold.  Drop distributions with weight less than T
- `U::Real` Merging threshold.  Merge if (μ_1-μ_2)^T * Σ^-1 * (μ_1-μ_2) < U
- `J_max::Integer` Maximum number of features
"""
function prune(x::GaussianMixture, T::Real, U::Real, J_max::Integer)
	l = 0
	J_k = x.N
	n = length(x.μ[1])

	# prune indices based on truncation threshold
	I = collect(1:J_k)
	I = I[x.w.>T]

	# setup output vectors
	w_new = zeros(J_k)
	μ_new = [zeros(n) for i=1:J_k]
	Σ_new = [zeros(n,n) for i=1:J_k]

	# perform pruning
	while length(I) != 0
		l = l+1
		j = I[argmax(x.w[I])]

		# Merge close features
		L = []
		for i in I
			delta = x.μ[i]-x.μ[j]
			if dot(delta, inv(x.Σ[i])*delta) < U
				push!(L,i)
			end
		end

		# determine merged parameters
		w_tilde = sum(x.w[L])
		μ_tilde = 1/w_tilde * sum(x.w[L].*x.μ[L])
		Σ_tilde = zeros(n,n)
		for i in L
			tmp = x.Σ[i] + (μ_tilde-x.μ[i])*(μ_tilde-x.μ[i])'
			Σ_tilde = Σ_tilde + x.w[i]*tmp
		end
		Σ_tilde = 1/w_tilde * Σ_tilde

		# store values
		w_new[l] = w_tilde
		μ_new[l] = μ_tilde
		Σ_new[l] = Σ_tilde

		# prune merged features
		I = setdiff(I,L)
	end


	if l > J_max 	# only keep the J_max with highest weights
		idxs = sort!([1:length(w_new);], by=i->(w_new[i]), rev=true)
		idxs = idxs[1:J_max]
		w_new = w_new[idxs]
		μ_new = μ_new[idxs]
		Σ_new = Σ_new[idxs]
	else		# drop zeros at end
		w_new = w_new[1:l]
		μ_new = μ_new[1:l]
		Σ_new = Σ_new[1:l]
	end

	x_new = GaussianMixture(w_new,μ_new,Σ_new)
	return x_new
end

### Utilities

"""
	MvNormalPDF(x::AbstractVector, μ::AbstractVector, Σ::AbstractMatrix)

Evaluate a multivariate normal pdf function parameterized by μ
and Σ at x

Arguments:
- `x::AbstractVector` Point to evaluate pdf
- `μ::AbstractVector` Mean of Multivariate Normal Gaussian
- `Σ::AbstractMatrix` Covariance Multivariate Normal Gaussian

Returns:
- `xp::AbstractVector` evaluation of multivariate normal pdf
"""
MvNormalPDF(x::AbstractVector, μ::AbstractVector, Σ::AbstractMatrix) =
	(2*pi)^(-length(μ)/2)*det(Σ)^(-1/2)*exp(-1/2*(x-μ)'*inv(Σ)*(x-μ))
