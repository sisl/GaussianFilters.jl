# Extended Kalman filter functions

"""
    predict(filter::ExtendedKalmanFilter, b0::GaussianBelief, u::AbstractVector)

Uses Extended Kalman filter to run prediction step on gaussian belief b0,
given control vector u.
"""
function predict(filter::ExtendedKalmanFilter, b0::GaussianBelief,
                u::AbstractVector{<:Number})

    # Motion update
    μp = predict(filter.d, b0.μ, u)
    F = jacobian(filter.d, b0.μ, u)

    Σp = F * b0.Σ * F' + filter.d.W
    return GaussianBelief(μp, Σp)
end

"""
    measure(filter::ExtendedKalmanFilter, bp::GaussianBelief, y::AbstractVector;
        u::AbstractVector = [false])

Uses Extended Kalman filter to run measurement update on predicted gaussian
belief bp, given measurement vector y. If u is specified and filter.o.D has
been declared, then matrix D will be factored into the y predictions.
"""
function measure(filter::ExtendedKalmanFilter, bp::GaussianBelief, y::AbstractVector{a};
                u::AbstractVector{b} = [false]) where {a<:Number, b<:Number}

    # Measurement update
    yp = measure(filter.o, bp.μ, u)
    H = jacobian(filter.o, bp.μ, u)

    # Kalman Gain
    K = bp.Σ * H' * inv(H * bp.Σ * H' + filter.o.V)

    # Measurement update
    μn = bp.μ + K * (y - yp)
    Σn = (I - K * H) * bp.Σ
    return GaussianBelief(μn, Σn)
end
