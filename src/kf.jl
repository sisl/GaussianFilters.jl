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
                filter::AbstractFilter;
                return_prediction = false) where {a<:Number, b<:Number}

    # predict
    bp = predict(b0, u, filter)

    # measure
    bn = measure(bp, y, filter; u = u)

    return return_prediction ? (bp, bn) : bn
end

# Kalman filter functions

"""
    predict(b0::GaussianBelief, u::Vector, filter::KalmanFilter)

Uses Kalman filter to run prediction step on gaussian belief b0, given control
vector u.
"""
function predict(b0::GaussianBelief, u::Vector{a},
            filter::KalmanFilter) where a<:Number

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

### Simulation functions ###

"""
    simulation(b0::GaussianBelief,action_sequence::Vector{Vector}},
        filter::AbstractFilter)

Run a simulation to get positions and measurements. Samples starting point from
GaussianBelief b0, the runs action_sequence with additive gaussian noise all
specified by AbstractFilter filter to return a simulated state and measurement
history.
"""
function simulation(b0::GaussianBelief, action_sequence::Vector{Vector{a}},
    filter::AbstractFilter) where a<:Number

    # make initial state
    s0 = b0.μ + cholesky(b0.Σ).L * randn(size(b0.Σ,1))

    # simulate action sequence
    state_history = [s0]
    measurement_history = Vector{Vector{typeof(s0[1])}}()
    for u in action_sequence
        xn, yn = simulate_step(state_history[end], u, filter)
        push!(state_history, xn)
        push!(measurement_history, yn)
    end

    ## TODO: Is it worth separating out start state from state_history
    return state_history, measurement_history
end

"""
    simulate_step(x::Vector, u::Vector, filter::KalmanFilter)

Run a step of simulation starting at state x, taking action u, and using the
motion and measurement equations specified by Kalman Filter filter.
"""
function simulate_step(x::Vector{a}, u::Vector{b},
    filter::KalmanFilter) where {a<:Number, b<:Number}

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
    run_filter(b0::GaussianBelief, action_history::Vector{Vector},
            measurement_history::Vector{Vector}, filter::AbstractFilter)

Given an initial belief b0, matched-size arrays for action and measurement
histories and a filter, update the beliefs using the filter, and return a
vector of all beliefs.
"""
function run_filter(b0::GaussianBelief, action_history::Vector{Vector{a}},
            measurement_history::Vector{Vector{b}},
            filter::AbstractFilter) where {a<:Number, b<:Number}

        # assert matching action and measurement sizes
        @assert length(action_history) == length(measurement_history)

        # initialize belief vector
        beliefs = [b0]

        # iterate through and update beliefs
        for (u, y) in zip(action_history, measurement_history)
            bn = update(beliefs[end], u, y, filter)
            push!(beliefs, bn)
        end

        return beliefs
end

"""
    beautify(belief_history::Vector{GaussianBelief};
        dims::Vector{Int}=[])

Given a history of beliefs, return a (time steps, state dim)-sized array of
predicted means and a (time steps, state dim, state dim)-sized array of
covariances. One can optionally specify dimensions indices dims to output
reduced state information.
"""
function beautify(belief_history::Vector{GaussianBelief{a,b}};
    dims::Vector{Int}=Vector{Int}()) where {a<:Number, b<:Number}

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
