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

### Simulation functions ###

"""
    simulation(filter::AbstractFilter, b0::GaussianBelief,
                action_sequence::Vector{AbstractVector}})

Run a simulation to get positions and measurements. Samples starting point from
GaussianBelief b0, the runs action_sequence with additive gaussian noise all
specified by AbstractFilter filter to return a simulated state and measurement
history.
"""
function simulation(filter::AbstractFilter, b0::GaussianBelief,
                    action_sequence::Vector{<:AbstractArray})

    # make initial state
    s0 = b0.μ + cholesky(b0.Σ).L * randn(size(b0.Σ,1))

    # simulate action sequence
    state_history = [s0]
    measurement_history = Vector{AbstractVector{typeof(s0[1])}}()
    for u in action_sequence
        xn, yn = simulate_step(filter, state_history[end], u)
        push!(state_history, xn)
        push!(measurement_history, yn)
    end

    ## TODO: Is it worth separating out start state from state_history
    return state_history, measurement_history
end

"""
    simulate_step(filter::KalmanFilter, x::AbstractVector, u::AbstractVector)

Run a step of simulation starting at state x, taking action u, and using the
motion and measurement equations specified by Kalman Filter filter.
"""
function simulate_step(filter::KalmanFilter, x::AbstractVector{<:Number},
                        u::AbstractVector{<:Number})

    # Linear Motion
    xn = filter.d.A * x + filter.d.B * u +
        cholesky(filter.d.W).L * randn(size(filter.d.W,1))

    # Linear Measurement
    yn = filter.o.C * xn +
        cholesky(filter.o.V).L * randn(size(filter.o.V,1))
    if !(filter.o.D[1,1] isa Bool)
        yn += filter.o.D * u
    end

    return xn, yn
end

"""
    run_filter(filter::AbstractFilter, b0::GaussianBelief, action_history::Vector{AbstractVector},
            measurement_history::Vector{AbstractVector})

Given an initial belief b0, matched-size arrays for action and measurement
histories and a filter, update the beliefs using the filter, and return a
vector of all beliefs.
"""
function run_filter(filter::AbstractFilter, b0::GaussianBelief, action_history::Vector{A},
            measurement_history::Vector{B}) where {A<:AbstractVector, B<:AbstractVector}

        # assert matching action and measurement sizes
        @assert length(action_history) == length(measurement_history)

        # initialize belief vector
        beliefs = [b0]

        # iterate through and update beliefs
        for (u, y) in zip(action_history, measurement_history)
            bn = update(filter, beliefs[end], u, y)
            push!(beliefs, bn)
        end

        return beliefs
end

"""
    unpack(belief_history::Vector{<:GaussianBelief};
        dims::Vector{Int}=[])

Given a history of beliefs, return an unpacked (time steps, state dim)-sized array of
predicted means and a (time steps, state dim, state dim)-sized array of
covariances. One can optionally specify dimensions indices dims to output
reduced state information.
"""
function unpack(belief_history::Vector{<:GaussianBelief};
    dims::Vector{Int}=Vector{Int}())

    # set default to condense all dimensions
    if length(dims) == 0
        dims = collect(1:length(belief_history[1].μ))
    end

    # set output sizes
    μ = zeros(typeof(belief_history[1].μ[1]),
        length(belief_history), length(dims))
    Σ = zeros(typeof(belief_history[1].μ[1]),
        length(belief_history), length(dims), length(dims))

    # iterate over belief_history and place elements appropriately
    for (i, belief) in enumerate(belief_history)
        μ[i,:] = belief.μ[dims]
        Σ[i,:,:] = belief.Σ[dims,dims]
    end

    return μ, Σ
end
### TODO: add in-place update!, predict!, and measure! functions
