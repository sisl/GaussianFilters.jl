Random.seed!(0)

# Dynamics
dt = 0.1
m = 50
A = [1 dt 0 0 ; 0 1 0 0 ; 0 0 1 dt; 0 0 0 1]
B = [0 0; dt/m 0; 0 0; 0 dt/m]
W = 0.1*Matrix{Float64}(I,4,4)

dmodel = LinearDynamicsModel(A,B,W);

# Measurement
C = [0 1.0 0 0; 0 0 0 1.0]
V = 0.5*Matrix{Float64}(I,2,2)

omodel = LinearObservationModel(C,V)

# Build KF
kf = KalmanFilter(dmodel,omodel)

# Simulation
b0 = GaussianBelief([0.0,0.0,0.0,0.0], 2.0*Matrix{Float64}(I,4,4))

times = 0:dt:10
Fmag = 1000
action_sequence = [[Fmag*cos(t), Fmag*sin(t)] for t in times]

sim_states, sim_measurements = simulation(kf, b0,action_sequence);

# Filter
filtered_beliefs = run_filter(kf, b0, action_sequence, sim_measurements)
μ, Σ = unpack(filtered_beliefs);

# Tests
@test kf.d isa LinearDynamicsModel
@test kf.o isa LinearObservationModel
@test length(sim_measurements) == length(action_sequence)

first_update = update(kf, b0, action_sequence[1], sim_measurements[1])
first_predict_measure = measure(kf, predict(kf, b0, action_sequence[1]),
    sim_measurements[1]; u = action_sequence[1])
array_isapprox(first_update.μ, first_predict_measure.μ)
array_isapprox(first_update.Σ, first_predict_measure.Σ)

@test length(filtered_beliefs[end].μ) == length(b0.μ)
@test size(filtered_beliefs[end].Σ) == size(b0.Σ)
@test size(μ) == (length(filtered_beliefs), length(filtered_beliefs[1].μ))
@test size(Σ) == (length(filtered_beliefs), size(filtered_beliefs[1].Σ)...)

# Test mathematical properties
@test all(isfinite.(sim_states[end]))  
@test all(isfinite.(sim_measurements[end])) 
@test length(sim_measurements[end]) == 2  
@test length(sim_states[end]) == 4 

# Test filtered belief properties
@test all(isfinite.(filtered_beliefs[end].μ))  
@test all(isfinite.(filtered_beliefs[end].Σ)) 
@test issymmetric(filtered_beliefs[end].Σ)  
@test isposdef(filtered_beliefs[end].Σ) 
@test tr(filtered_beliefs[end].Σ) < 30.0 
