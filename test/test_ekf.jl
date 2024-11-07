Random.seed!(0)

# Dynamics
dt = 0.001
J = [1; 5; 5]
function dynamics(x,u)
    xp = x + dt*u./J
    xp[1] += dt*(J[2]-J[3])/J[1]*x[2]*x[3]
    xp[2] += dt*(J[3]-J[1])/J[2]*x[3]*x[1]
    xp[3] += dt*(J[1]-J[2])/J[3]*x[1]*x[2]
    return xp
end

W = 0.001*Matrix{Float64}(I,3,3)

dmodel = NonlinearDynamicsModel(dynamics,W);

# Measurement
c=10
function observe(x,u)
    y = x
    y = min.(y,c)
    y = max.(y,-c)
    return y
end

V = 0.3*Matrix{Float64}(I,3,3)

omodel = NonlinearObservationModel(observe,V)

# Build EKF
ekf = ExtendedKalmanFilter(dmodel,omodel)

# Simulation
times = 0:dt:0.1
action_sequence = [[0.0,0.0,0.0] for t in times]

b0 = GaussianBelief([10.0,0.0,0.0], Matrix{Float64}(I,3,3))

sim_states, sim_measurements = simulation(ekf,b0,action_sequence)

# Filter
filtered_beliefs = run_filter(ekf, b0, action_sequence, sim_measurements)
μ, Σ = unpack(filtered_beliefs)

# Tests
@test ekf.d isa NonlinearDynamicsModel
@test ekf.o isa NonlinearObservationModel
@test length(sim_measurements) == length(action_sequence)


first_update = update(ekf, b0, action_sequence[1], sim_measurements[1])
first_predict = predict(ekf, b0, action_sequence[1])
first_predict_measure = measure(ekf, first_predict, sim_measurements[1], u = action_sequence[1])
array_isapprox(first_update.μ, first_predict_measure.μ)
array_isapprox(first_update.Σ, first_predict_measure.Σ)


@test length(filtered_beliefs[end].μ) == length(b0.μ)
@test size(filtered_beliefs[end].Σ) == size(b0.Σ)
@test size(μ) == (length(filtered_beliefs), length(filtered_beliefs[1].μ))
@test size(Σ) == (length(filtered_beliefs), size(filtered_beliefs[1].Σ)...)

# Test mathematical properties
@test all(isfinite.(sim_states[end]))  
@test all(isfinite.(sim_measurements[end])) 
@test length(sim_measurements[end]) == 3  
@test length(sim_states[end]) == 3 

# Test filtered belief properties
@test all(isfinite.(filtered_beliefs[end].μ))  
@test all(isfinite.(filtered_beliefs[end].Σ)) 
@test issymmetric(filtered_beliefs[end].Σ)  
@test isposdef(filtered_beliefs[end].Σ) 
@test tr(filtered_beliefs[end].Σ) < 1.0 
