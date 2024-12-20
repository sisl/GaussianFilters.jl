# Generic update function

"""
    update(filter::AbstractFilter, b0::GaussianBelief, u::AbstractVector,
        y::AbstractVector)

Uses AbstractFilter filter to update gaussian belief b0, given control vector
u and measurement vector y.
"""
function update(filter::AbstractFilter, b0::GaussianBelief,
                u::AbstractVector{<:Number}, y::AbstractVector{<:Number})

    # predict
    bp = predict(filter, b0, u)

    # measure
    bn = measure(filter, bp, y; u = u)

    return bn
end

function update_with_info(filter::AbstractFilter, b0::GaussianBelief,
    u::AbstractVector{<:Number}, y::AbstractVector{<:Number})

    # predict
    bp = predict(filter, b0, u)

    # measure
    bn, info = measure_info(filter, bp, y; u = u)

    return bn, info
end

# Kalman filter functions

"""
    predict(filter::KalmanFilter, b0::GaussianBelief, u::AbstractVector)

Uses Kalman filter to run prediction step on gaussian belief b0, given control
vector u.
"""
function predict(filter::KalmanFilter, b0::GaussianBelief,
            u::AbstractVector{<:Number})

    # Motion update
    μp = filter.d.A * b0.μ + filter.d.B * u
    Σp = filter.d.A * b0.Σ * filter.d.A' + filter.d.W
    return GaussianBelief(μp, Σp)
end

"""
    measure(filter::KalmanFilter, bp::GaussianBelief, y::AbstractVector;
        u::AbstractVector = [false])

Uses Kalman filter to run measurement update on predicted gaussian belief bp,
given measurement vector y. If u is specified and filter.o.D has been declared,
then matrix D will be factored into the y predictions
"""
function measure(filter::KalmanFilter, bp::GaussianBelief, y::AbstractVector{<:Number};
                u::AbstractVector{<:Number} = [false])
    return first(measure_info(filter, bp, y; u = u))
end

function measure_info(filter::KalmanFilter, bp::GaussianBelief, y::AbstractVector{<:Number};
    u::AbstractVector{<:Number} = [false])

    # Kalman Gain
    ΣY = filter.o.C * bp.Σ * filter.o.C' + filter.o.V
    K = bp.Σ * filter.o.C' * inv(ΣY)

    # Predicted measurement
    yp = measure(filter.o, bp.μ, u)

    # Measurement update
    μn = bp.μ + K * (y-yp)
    Σn = (I - K * filter.o.C) * bp.Σ
    
    info = (innovation_cov = ΣY, kalman_gain = K, predicted_measurement = yp)

    return GaussianBelief(μn, Σn), info
end
