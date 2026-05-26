using GaussianFilters
using POMDPs
using LinearAlgebra
using Random

# This example shows that a GaussianFilters.AbstractFilter doubles as a
# POMDPs.jl belief updater. Once POMDPs is loaded, the package extension
# wires up POMDPs.update(filter, b, a, o) -> GaussianFilters.update(...).

# A simple linear system: position-velocity 1D motion observed in position.
dt = 0.1
A = [1.0 dt; 0.0 1.0]
B = reshape([0.0, dt], 2, 1)
W = 0.01 * Matrix{Float64}(I, 2, 2)
dmodel = LinearDynamicsModel(A, B, W)

C = [1.0 0.0]
V = 0.05 * Matrix{Float64}(I, 1, 1)
omodel = LinearObservationModel(C, V)

kf = KalmanFilter(dmodel, omodel)
b0 = GaussianBelief([0.0, 0.0], Matrix{Float64}(I, 2, 2))

# Simulate a few steps so we have an action/observation history.
rng = MersenneTwister(0)
sim_states, sim_observations = simulation(kf, b0, [[1.0] for _ in 1:5], rng)

# Use POMDPs.update (which dispatches to GaussianFilters.update through
# the package extension) to filter the observations.
belief = POMDPs.initialize_belief(kf, b0)
for (a, o) in zip([[1.0] for _ in 1:5], sim_observations)
    global belief = POMDPs.update(kf, belief, a, o)
end

println("Final belief mean: ", belief.μ)
println("Final belief covariance: ")
display(belief.Σ)
