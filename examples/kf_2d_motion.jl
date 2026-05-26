using GaussianFilters

using LinearAlgebra
using Distributions
using Random

Random.seed!(1);

dt = 0.1
m = 50
A = [1 dt 0 0 ; 0 1 0 0 ; 0 0 1 dt; 0 0 0 1]
B = [0 0; dt/m 0; 0 0; 0 dt/m]
W = 0.1*Matrix{Float64}(I,4,4)

# build linear dynamics model
dmodel = LinearDynamicsModel(A,B,W);

C = [0 1.0 0 0; 0 0 0 1.0]
V = 0.5*Matrix{Float64}(I,2,2)

# build linear observation model
omodel = LinearObservationModel(C,V)

# build kf
kf = KalmanFilter(dmodel,omodel);

b0 = GaussianBelief([0.0,0.0,0.0,0.0], 2.0*Matrix{Float64}(I,4,4))

times = 0:dt:10
Fmag = 1000
action_sequence = [[Fmag*cos(t), Fmag*sin(t)] for t in times]

sim_states, sim_measurements = simulation(kf,b0,action_sequence);

filtered_beliefs = run_filter(kf, b0, action_sequence, sim_measurements)

# turn array of belief structs into simple tensors.
μ, Σ = unpack(filtered_beliefs;dims=[1,3]);

using Plots

x = [x[1] for x in sim_states]
y = [x[3] for x in sim_states]
plot(x, y)

plot!(μ[:,1], μ[:,2])

for i in 1:20:size(μ,1)
    x,y = belief_ellipse(μ[i,:], Σ[i,:,:])
    plot!(x,y)
end
plot!(legend=false)

vx = [x[2] for x in sim_states]
vy = [x[4] for x in sim_states]
plot(vx, vy)
