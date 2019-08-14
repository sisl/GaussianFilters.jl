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
function unscented_transform(b::GaussianBelief, λ::Int=2, α::Number=1,
    β::Number=0; decomp_method::String = "cholesky")

    # compute state space dimensionality
    n = length(b.μ)

    # compute weights
    w_μ = 1/(2.0*(n+λ))*ones(2*n+1)
    w_μ[1] = λ/(λ+n)
    w_Σ = copy(w_m)
    w_Σ[1] += (1 - α^2 + β) # Per ProbRob formulation

    # compute Σ^0.5 using different methods
    if decomp_method == "cholesky" # default
        s = cholesky(b.Σ).L
    else if decomp_method == "unit-axes"
        # TODO: SVD decomp method
        error("Unit axes method not yet implemented.")
    else if decomp_method == "tilt-axes"
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
    unscented_transform_inverse(points::Array{Array}, w_μ::Array, w_Σ::Array)

Convert a 2n+1 sigma points and weights back to a single measure for
mean and covariance (GaussianBelief).

Uses formulation from ProbRob for α/β parameters on a separate covariance
weighting, although this is not necessary.

"""
function unscented_transform_inverse(points::Array{Array{a}}, w_m::Array{b},
    w_c::Array{c}) where {a<:Number, b<:Number, c<:Number}

    # calculate weighted mean
    μ = sum(points .* w_μ)

    # calculated weighted covariance
    diff = hcat(points...) .- μ
    scaled = diff .* sqrt.(w_c)'
    Σ = scaled * scaled'

    return GaussianBelief(μ, Σ)
end

#=
"""
Unscented Kalman Update

Rework into separate predict and measure functions, with update=predict+measure
"""
function UnscentedKalmanUpdate(μ, Σ, u, y, g, h, W, V ;dt=0.001)
    n = length(μ)

    # predict
    points, w_m, w_c = UT(μ,Σ)
    pointsp = zeros(size(points))
    for i in 1:(2*n+1)
        pointsp[:,i] = g(points[:,i], u;dt=dt)
    end
    μ_pred, Σb_pred = invUT(pointsp, w_m, w_c)
    Σ_pred = Σb_pred + W

    # update
    points, w_m, w_c = UT(μ_pred,Σ_pred)
    y_preds = zeros(length(y),(2*n+1))
    for i in 1:(2*n+1)
        y_preds[:,i] = h(points[:,i])
    end
    y_pred = y_preds*w_m

    scaled_ydiff = (y_preds .- y_pred).*sqrt.(w_c)'
    Σ_Y = scaled_ydiff*scaled_ydiff' + V

    scaled_xdiff = (points .- μ_pred).*sqrt.(w_c)'
    Σ_XY = scaled_xdiff*scaled_ydiff'

    μ = μ_pred + Σ_XY*inv(Σ_Y)*(y-y_pred)
    Σ = Σ_pred - Σ_XY*inv(Σ_Y)*Σ_XY'
    Σ = Symmetric(Σ)
    return μ, Σ
end
=#

# Unscented Kalman Filter functions

"""
    predict(b0::GaussianBelief, u::Vector, filter::UnscentedKalmanFilter)

Uses Extended Kalman filter to run prediction step on gaussian belief b0,
given control vector u.
"""
function predict(b0::GaussianBelief, u::Vector{a},
            filter::UnscentedKalmanFilter) where a<:Number

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
    measure(bp::GaussianBelief, y::Vector, filter::UnscentedKalmanFilter;
        u::Vector = [false])

Uses Extended Kalman filter to run measurement update on predicted gaussian
belief bp, given measurement vector y. If u is specified and filter.o.D has
been declared, then matrix D will be factored into the y predictions.
"""
function measure(bp::GaussianBelief, y::Vector{a}, filter::UnscentedKalmanFilter;
                u::Vector{b} = [false]) where {a<:Number, b<:Number}

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
    simulate_step(x::Vector, u::Vector, filter::UnscentedKalmanFilter)

Run a step of simulation starting at state x, taking action u, and using the
motion and measurement equations specified by Kalman Filter filter.
"""
function simulate_step(x::Vector{a}, u::Vector{b},
    filter::UnscentedKalmanFilter) where {a<:Number, b<:Number}

    # Motion

    xn = cholesky(filter.d.W).L * randn(size(filter.d.W,1))

    # Linear Motion
    if filter.d isa LinearDynamicsModel
        xn += filter.d.A * x + filter.d.B * u

    # Nonlinear Motion
    elseif filter.d isa NonlinearDynamicsModel
        xn += filter.d.f(x,u)

    else
        error("Unsupported UKF Dynamics Model Type: " *
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
        error("Unsupported UKF Observation Model Type: " *
            string(typeof(filter.o)))
    end

    return xn, yn
end
