using GaussianFilters

using LinearAlgebra
using Distributions
using Random

Random.seed!(2019);

# nonlinear dynamics function. must be a function of both states (x) and actions (u) even if either are not used.
dt=0.1
function step(x,u)
    v,ϕ = u
    xp = x
    xp += dt*([v*cos(x[3]);v*sin(x[3]); ϕ])
    xp[3] = xp[3]
    return xp
end  

W = dt*0.05*Matrix{Float64}(I,3,3)

# build dynamics model
dmodel = NonlinearDynamicsModel(step,W);

# nonlinear observation function. must be a function of both states (x) and actions (u) even if either are not used.
function observe(x,u)
    y = [norm(x[1:2])]
end

V = 0.1*Matrix{Float64}(I,1,1)

# build observation model
omodel = NonlinearObservationModel(observe,V)

# build ukf
ukf = UnscentedKalmanFilter(dmodel,omodel);

times = 0:dt:40
action_sequence = [[1.0, 2*sin(pi/2*t)] for t in times]

b0 = GaussianBelief([0.0,0.0,0.0], 0.01*Matrix{Float64}(I,3,3))

sim_states, sim_measurements = simulation(ukf,b0,action_sequence);

filtered_beliefs = run_filter(ukf, b0, action_sequence, sim_measurements)

# turn array of belief structs into simple tensors.
μ, Σ = unpack(filtered_beliefs);

using Plots

plot(times,hcat(sim_states[2:end]...)',labels=["px","py","theta"])

plot!(times,μ[2:end,:], linestyle=:dash)
plot!(legend=false)

# form and run ekf
ekf = ExtendedKalmanFilter(dmodel,omodel);
filtered_beliefs_ekf = run_filter(ekf, b0, action_sequence, sim_measurements)

# turn array of belief structs into simple tensors.
μ_ekf, Σ_ekf = unpack(filtered_beliefs_ekf); 

plot!(times,μ_ekf[2:end,:], linestyle=:dot)
plot!(legend=false)
