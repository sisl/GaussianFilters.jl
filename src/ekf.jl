"""
    Work into extended Kalman update, separate by predict, measure,
    update where update = predict, measure, use autodiff to compute G and H
"""

function ExtendedKalmanUpdate(μ, Σ, u, y, g, h, W, V;dt=0.001)
    μ_pred = g(μ, u;dt=dt)
    # G = dg/dmu
    Σ_pred = G*Σ*G' + W

    # H = dh/dmu
    K = Σ_pred*H'*inv(H*Σ_pred*H'+V)
    μ = μ_pred + K*(y-h(μ_pred))
    Σ = (I - K*H)*Σ_pred
    return μ, Σ
end


# Extended Kalman filter functions

"""
    predict(b0::GaussianBelief, u::Vector, filter::ExtendedKalmanFilter)

Uses Extended Kalman filter to run prediction step on gaussian belief b0,
given control vector u.
"""
function predict(b0::GaussianBelief, u::Vector{a},
            filter::ExtendedKalmanFilter) where a<:Number

    # Motion update
    μp = filter.d.A * b0.μ + filter.d.B * u
    Σp = filter.d.A * b0.Σ * filter.d.A' + filter.d.W
    return GaussianBelief(μp, Σp)
end

"""
    measure(bp::GaussianBelief, y::Vector, filter::ExtendedKalmanFilter;
        u::Vector = [false])

Uses Extended Kalman filter to run measurement update on predicted gaussian
belief bp, given measurement vector y. If u is specified and filter.o.D has
been declared, then matrix D will be factored into the y predictions.
"""
function measure(bp::GaussianBelief, y::Vector{a}, filter::ExtendedKalmanFilter;
                u::Vector{b} = [false]) where {a<:Number, b<:Number}
    # Kalman Gain
    K = bp.Σ * filter.o.C' *
        inv(filter.o.C * bp.Σ * filter.o.C' + filter.o.V)

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

### Simulation function ###
"""
    simulate_step(x::Vector, u::Vector, filter::ExtendedKalmanFilter)

Run a step of simulation starting at state x, taking action u, and using the
motion and measurement equations specified by Kalman Filter filter.
"""
function simulate_step(x::Vector{a}, u::Vector{b},
    filter::ExtendedKalmanFilter) where {a<:Number, b<:Number}

    # Motion

    xn = cholesky(filter.d.W).L * randn(size(filter.d.W,1))

    # Linear Motion
    if filter.d isa LinearDynamicsModel
        xn += filter.d.A * x + filter.d.B * u

    # Nonlinear Motion
    elseif filter.d isa NonlinearDynamicsModel
        xn += filter.d.f(x,u)

    else
        error("Unsupported EKF Dynamics Model Type: " *
            string(typeof(filter.d)))
    end

    # Measurement

    yn = cholesky(filter.o.V).L * randn(size(filter.o.V,1))

    # Linear Measurement
    if filter.o isa LinearObservationModel
        yn += filter.o.C * xn
        if !(filter.o.D[1,1] isa Bool)
            yn += filter.o.D * u
        end

    # Nonlinear Measurement
    elseif filter.o isa NonlinearObservationModel
        yn += filter.o.h(xn,u)
    else
        error("Unsupported EKF Observation Model Type: " *
            string(typeof(filter.o)))
    end

    return xn, yn
end
