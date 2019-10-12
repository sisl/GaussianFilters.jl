### SIMULATION FUNCTIONS ###

"""
    simulate_step(filter::AbstractFilter, x::AbstractVector, u::AbstractVector, rng::AbstractRNG=Random.GLOBAL_RNG)

Run a step of simulation starting at state x, taking action u, and using the
motion and measurement equations specified by the filter.
"""
function simulate_step(filter::AbstractFilter, x::AbstractVector{<:Number},
                        u::AbstractVector{<:Number}, rng::AbstractRNG=Random.GLOBAL_RNG)
    
    xn = predict(filter.d, x, u, rng)
    
    yn = measure(filter.o, xn, u, rng)

    return xn, yn
end


"""
    simulation(filter::AbstractFilter, b0::GaussianBelief,
                action_sequence::Vector{AbstractVector}})

Run a simulation to get positions and measurements. Samples starting point from
GaussianBelief b0, the runs action_sequence with additive gaussian noise all
specified by AbstractFilter filter to return a simulated state and measurement
history.
"""
function simulation(filter::AbstractFilter, b0::GaussianBelief,
                    action_sequence::Vector{<:AbstractArray},
                    rng::AbstractRNG = Random.GLOBAL_RNG)

    # make initial state
    s0 = rand(rng, b0)

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
