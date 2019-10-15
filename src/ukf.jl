# Sigma point transform functions

"""
unscented_transform(b::GaussianBelief, λ::Int=2, α::Number=1,
    β::Number=0; decomp_method::String = "cholesky")

Convert a single GaussianBelief (mean and covariance) to set of 2n+1 sigma
points, with n being the dimensionality of the state space. Return an array of
the points, and arrays for weights used in mean and covariance calculations.
Uses formulation from ProbRob for α/β parameters on a separate covariance
weighting, although this is not necessary.

"""
function unscented_transform(b::GaussianBelief, λ::Number=2, α::Number=1,
    β::Number=0; decomp_method::String = "cholesky")

    # compute state space dimensionality
    n = length(b.μ)

    # compute weights
    w_μ = 1/(2.0*(n+λ))*ones(2*n+1)
    w_μ[1] = λ/(λ+n)
    w_Σ = copy(w_μ)
    w_Σ[1] += (1 - α^2 + β) # Per ProbRob formulation

    # compute Σ^0.5 using different methods
    if decomp_method == "cholesky" # default
        s = cholesky(b.Σ).L
    elseif decomp_method == "unit-axes"
        # TODO: SVD decomp method
        error("Unit axes method not yet implemented.")
    elseif decomp_method == "tilt-axes"
        # TODO: SVD decomp method with tilt transformation
        error("Tilt axes method not yet implemented.")
    else
        error("Undefined unscented transform Σ sqrt.")
    end

    # compute sigma points
    points = []
    push!(points,b.μ)
    for i in 1:n
        push!(points, b.μ + sqrt(n+λ)*s[:,i])
        push!(points, b.μ - sqrt(n+λ)*s[:,i])
    end

    return points, w_μ, w_Σ
end

"""
    unscented_transform_inverse(points::Vector{AbstractVector}, w_μ::Vector,
        w_Σ::Vector)

Convert a 2n+1 sigma points and weights back to a single measure for
mean and covariance (GaussianBelief).

Uses formulation from ProbRob for α/β parameters on a separate covariance
weighting, although this is not necessary.

"""
function unscented_transform_inverse(points::Vector{A}, w_μ::Vector{B},
    w_Σ::Vector{C}) where {A<:AbstractVector, B<:Number, C<:Number}

    # calculate weighted mean
    μ = sum(points .* w_μ)

    # calculated weighted covariance
    diff = hcat(points...) .- μ
    scaled = diff .* sqrt.(w_Σ)'
    Σ = scaled * scaled'

    return GaussianBelief(μ, Σ)
end

# Unscented Kalman Filter functions

"""
    predict(filter::UnscentedKalmanFilter, b0::GaussianBelief, u::AbstractVector)

Uses Unscented Kalman filter to run prediction step on gaussian belief b0,
given control vector u.
"""
function predict(filter::UnscentedKalmanFilter, b0::GaussianBelief, u::AbstractVector{<:Number})

    # Motion update

    n = length(b0.μ)

    # approximate Gaussian belief with sigma points
    points, w_μ, w_Σ = unscented_transform(b0, filter.λ, filter.α, filter.β)

    # iterate over each sigma point and propagate it through motion function
    pointsp = [predict(filter.d, point, u) for point in points]

    # apply inverse unscented transform to approximate new Gaussian
    bp = unscented_transform_inverse(pointsp, w_μ, w_Σ)

    # add process noise
    Σp = bp.Σ + filter.d.W

    return GaussianBelief(bp.μ, Σp)
end

"""
    measure(filter::UnscentedKalmanFilter, bp::GaussianBelief, y::AbstractVector;
        u::AbstractVector = [false])

Uses Unscented Kalman filter to run measurement update on predicted gaussian
belief bp, given measurement vector y. If u is specified and filter.o.D has
been declared, then matrix D will be factored into the y predictions.
"""
function measure(filter::UnscentedKalmanFilter, bp::GaussianBelief, y::AbstractVector{a};
                u::AbstractVector{b} = [false]) where {a<:Number, b<:Number}

    # Measurement update

    # approximate Gaussian belief with sigma points
    points, w_μ, w_Σ = unscented_transform(bp, filter.λ, filter.α, filter.β)

    # iterate over sigma points and computed expected measurement
    ysp = [measure(filter.o, point, u) for point in points]
    
    # compute expected measurement
    yp = sum(ysp .* w_μ)

    # compute marginal covariance components
    ydiff = hcat(ysp...) .- yp
    scaled_ydiff = ydiff .* sqrt.(w_Σ)'
    Σ_Y = scaled_ydiff * scaled_ydiff' + filter.o.V

    xdiff = hcat(points...) .- bp.μ
    scaled_xdiff = xdiff .* sqrt.(w_Σ)'
    Σ_XY = scaled_xdiff * scaled_ydiff'

    # measurement update
    μn = bp.μ + Σ_XY * inv(Σ_Y) * (y - yp)
    Σn = bp.Σ - Σ_XY * inv(Σ_Y) * Σ_XY'
    return GaussianBelief(μn, Σn)
end
