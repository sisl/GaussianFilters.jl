# Extended Kalman filter functions

"""
    predict(filter::ExtendedKalmanFilter, b0::GaussianBelief, u::AbstractVector)

Uses Extended Kalman filter to run prediction step on gaussian belief b0,
given control vector u.
"""
function predict(filter::ExtendedKalmanFilter, b0::GaussianBelief,
                u::AbstractVector{<:Number})

    # Motion update

    # Linear motion
    if filter.d isa LinearDynamicsModel
        μp = filter.d.A * b0.μ + filter.d.B * u
        F = filter.d.A

    # Nonlinear motion
    elseif filter.d isa NonlinearDynamicsModel

        # Nonlinear predicted motion
        μp = filter.d.f(b0.μ, u)

        # Nonlinear Jacobian calculated with ForwardDiff.jl
        F = ForwardDiff.jacobian(μ -> filter.d.f(μ, u), b0.μ)

    else
        error("Unsupported EKF Dynamics Model Type: " *
            string(typeof(filter.d)))
    end

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

    # Linear measurement
    if filter.o isa LinearObservationModel

        # Linear predicted measurement
        yp = filter.o.C * bp.μ
        if !(filter.o.D[1,1] isa Bool)
            if u[1]==false
                @warn "D matrix specified in measurement model but not being used"
            else
                yp = yp + filter.o.D * u
            end
        end

        # Linear Jacobian
        H = filter.o.C

    # Nonlinear measurement
    elseif filter.o isa NonlinearObservationModel

        # Nonlinear predicted measurement
        yp = filter.o.h(bp.μ, u)

        # Nonlinear Jacobian calculated with ForwardDiff.jl
        H = ForwardDiff.jacobian(μ -> filter.o.h(μ, u), bp.μ)

    else
        error("Unsupported EKF Observation Model Type: " *
            string(typeof(filter.o)))
    end

    # Kalman Gain
    K = bp.Σ * H' * inv(H * bp.Σ * H' + filter.o.V)

    # Measurement update
    μn = bp.μ + K * (y - yp)
    Σn = (I - K * H) * bp.Σ
    return GaussianBelief(μn, Σn)
end

### Simulation function ###
"""
    simulate_step(filter::ExtendedKalmanFilter, x::AbstractVector, u::AbstractVector)

Run a step of simulation starting at state x, taking action u, and using the
motion and measurement equations specified by Kalman Filter filter.
"""
function simulate_step(filter::ExtendedKalmanFilter, x::AbstractVector{<:Number},
    u::AbstractVector{<:Number})

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
