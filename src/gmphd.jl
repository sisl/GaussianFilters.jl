"""
	update(phd::PHDFilter, b0::GaussianMixture, Z::Vector{AbstractVector},
		T::Real, U::Real, J_max::Integer)

Perform an update step using a GMPHD Filter.

Arguments:
- `phd::PHDFilter` PHD filter to step through. [PHD Filter]
- `b0::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Matrix` Matrix of measurements. [Measurements]
- `T::Real` Truncation threshold.  Drop distributions with weight less than T
- `U::Real` Merging threshold.  Merge if (μ_1-μ_2)^T * Σ^-1 * (μ_1-μ_2) < U
- `J_max::Integer` Maximum number of features

Returns:
- `bn::GaussianMixture` Describe first return value. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function update(phd::PHDFilter, b0::GaussianMixture,
	Z::Vector{<:AbstractVector{<:Number}}, T::Real, U::Real, J_max::Integer)

	bp = predict(phd, b0)
	bm = measure(phd, bp, Z)
    bn = prune(bm, T, U, J_max)
    return bn
end

"""
	predict(phd::PHDFilter, b0::GaussianMixture)

Make a prediction on next state based on PHD dynamics

Arguments:
- `phd::PHDFilter` PHD filter to step through. [PHD Filter]
- `b0::GaussianMixture` Prior state distribution. [Gaussian Mixture]

Returns:
- `bp::GaussianMixture` Predicted next state distribution. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function predict(phd::PHDFilter, b0::GaussianMixture)

    J_k = b0.N
    n = length(b0.μ[1])
    J_γ = phd.γ.N
    J_β = phd.spawn.β.N
    J_dyn = length(phd.dyn)
    tot = J_γ + J_β*J_k + J_k*J_dyn
    w_p = [0.0 for i in 1:tot]
    μ_p = [zeros(n) for i in 1:tot]
    Σ_p = [zeros(n,n) for i in 1:tot]

    # 1. Prediction for birth targets
    i = 0
    for j in 1:J_γ
        i += 1
        w_p[i] = phd.γ.w[j]
        μ_p[i] = phd.γ.μ[j]
        Σ_p[i] = phd.γ.Σ[j]
    end


    for j in 1:J_β, ℓ in 1:J_k
        i += 1
        spawn_dyn = phd.spawn.dyn[j]
        w_p[i] = b0.w[ℓ]*phd.spawn.β.w[j]
        μ_p[i] =  spawn_dyn.d + spawn_dyn.A*b0.μ[ℓ]
        Σ_p[i] = spawn_dyn.Q + spawn_dyn.A*b0.Σ[ℓ]*spawn_dyn.A'
    end

    # 2. Prediction for existing targets
    for j in 1:J_k
        for k in 1:J_dyn
            i += 1
            w_p[i] = phd.Ps*b0.w[j]
            μ_p[i] = phd.dyn[k].A*b0.μ[j]
            Σ_p[i] = phd.dyn[k].Q + phd.dyn[k].A*b0.Σ[j]*phd.dyn[k].A'
        end
    end
    @assert i == tot "Poor Indexing"
    bp = GaussianMixture(w_p, μ_p, Σ_p)
	return bp
end

"""
	measure(phd::PHDFilter, bp::GaussianMixture, Z::Vector{AbstractVector})

Perform a measurement update on a predicted PHD next state.

Arguments:
- `phd::PHDFilter` PHD filter to step through. [PHD Filter]
- `bp::GaussianMixture` Prior state distribution. [Gaussian Mixture]
- `Z::Vector{AbstractVector}` Array of measurements. [Measurements]

Returns:
- `bm::GaussianMixture` Describe first return value. [Gaussian Mixture]

References:
1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis
density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.
"""
function measure(phd::PHDFilter, bp::GaussianMixture,
				Z::Vector{<:AbstractVector{<: Real}})

	# 3. Construction of PHD update components

    J_p = bp.N

    η_n = []
    S_n = []
    K_n = []
    P_n = []
    for j in 1:J_p
        push!(η_n, phd.meas.C*bp.μ[j])

        S_add = phd.meas.R + phd.meas.C*bp.Σ[j]*phd.meas.C'
        push!(S_n, S_add)

        K_add = bp.Σ[j]*phd.meas.C'*inv(S_add)
        push!(K_n, K_add)

        P_add = (1.0I - K_add*phd.meas.C)*bp.Σ[j]
        push!(P_n, P_add)
    end

    # 4. Update
    w_n = similar(bp.w,0)
    μ_n = similar(bp.μ,0)
    Σ_n = similar(bp.Σ,0)
    for j in 1:J_p
        push!(w_n,(1-phd.Pd)*bp.w[j])
        push!(μ_n, bp.μ[j])
        push!(Σ_n, bp.Σ[j])
    end

    for z in Z
        w_ntemp = []
        for j in 1:J_p
            push!(w_ntemp, phd.Pd*bp.w[j]*MvNormalPDF(z,η_n[j],S_n[j]))
            push!(μ_n, bp.μ[j] + K_n[j]*(z-η_n[j]))
            push!(Σ_n, P_n[j])
        end
        append!(w_n, w_ntemp ./ (phd.κ(z) + sum(w_ntemp)))
    end

    xm = GaussianMixture(w_n, μ_n, Σ_n)
    return xm
end

"""
    prune(b::GaussianMixture, T::Real, U::Real, J_max::Integer)

Prune the posterior Gaussian-Mixture next PHD state

Arguments:
- `b::GaussianMixture` Posterior state distribution. [Gaussian Mixture]
- `T::Real` Truncation threshold.  Drop distributions with weight less than T
- `U::Real` Merging threshold.  Merge if (μ_1-μ_2)^T * Σ^-1 * (μ_1-μ_2) < U
- `J_max::Integer` Maximum number of features
"""
function prune(b::GaussianMixture, T::Real, U::Real, J_max::Integer)
	l = 0
	J_k = b.N
	n = length(b.μ[1])

	# prune indices based on truncation threshold
	I = collect(1:J_k)
	I = I[b.w.>T]

	# setup output vectors
	w_new = zeros(J_k)
	μ_new = [zeros(n) for i=1:J_k]
	Σ_new = [zeros(n,n) for i=1:J_k]

	# perform pruning
	while length(I) != 0
		l = l+1
		j = I[argmax(b.w[I])]

		# Merge close features
		L = []
		for i in I
			delta = b.μ[i]-b.μ[j]
			if dot(delta, inv(b.Σ[i])*delta) < U
				push!(L,i)
			end
		end

		# determine merged parameters
		w_tilde = sum(b.w[L])
		μ_tilde = 1/w_tilde * sum(b.w[L].*b.μ[L])
		Σ_tilde = zeros(n,n)
		for i in L
			tmp = b.Σ[i] + (μ_tilde-b.μ[i])*(μ_tilde-b.μ[i])'
			Σ_tilde = Σ_tilde + b.w[i]*tmp
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

	bn = GaussianMixture(w_new,μ_new,Σ_new)
	return bn
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
