# Generic update function

"""
update(b0::GaussianBelief, u::Vector, y::Vector, filter::AbstractFilter;
    return_prediction = false)

Uses AbstractFilter filter to update gaussian belief b0, given control vector
u and measurement vector y. If return_preduction is set to true, update
also returns the predicted state (before the measurement update) as a second
output
"""
function update(b0::GaussianBelief, u::Vector{a}, y::Vector{b},
                filter::AbstractFilter; return_prediction = false)
                where {a<:Number, b<:Number}

    # predict
    bp = predict(b0, u, filter)

    # measure
    bn = measure(bp, y, filter; u = u)

    return_prediction ? return bp, bn : return bn
end

# Kalman filter functions

"""
predict(b0::GaussianBelief, u::Vector, filter::KalmanFilter)

Uses Kalman filter to run prediction step on gaussian belief b0, given control
vector u.
"""
function predict(b0::GaussianBelief, u::Vector{a}, filter::KalmanFilter)
                where a<:Number

    # Motion update
    μp = filter.d.A * b0.μ + filter.d.B * u
    Σp = filter.d.A * b0.Σ * filter.d.A' + filter.d.W
    return GaussianBelief(μp, Σp)
end

"""
measure(bp::GaussianBelief, y::Vector, filter::KalmanFilter;
    u::Vector = [false])

Uses Kalman filter to run measurement update on predicted gaussian belief bp,
given measurement vector y. If u is specified and filter.o.D has been declared,
then matrix D will be factored into the y predictions
"""
function measure(bp::GaussianBelief, y::Vector{a}, filter::KalmanFilter;
                u::Vector{b} = [false]) where a<:Number
    # Kalman Gain
    K = bp.Σ * filter.o.C'
        * inv(filter.o.C * bp.Σ * filter.o.C'+V)

    # Predicted measurement
    yp = filter.o.C * bp.μ
    if !(filter.o.D[1,1] isa Bool)
        if u[1]==false
            @warn "D matrix specified in measurement model but not being used"
        else
            yp = yp + filter.o.D * u
        end
    end

    # Measurement update
    μn = bp.μ + K * (y-yp)
    Σn = (I - K * filter.o.C) * bp.Σ
    return GaussianBelief(μn, Σn)
end

### TODO: add in place update!, predict!, and measure! functions
