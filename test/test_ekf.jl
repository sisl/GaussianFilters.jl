rng = StableRNG(0)

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

sim_states, sim_measurements = simulation(ekf,b0,action_sequence, rng)

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

array_isapprox(sim_states[end], [9.476812080278332, -0.31100645773037805, 0.6082138148099138]; atol=0.001)
array_isapprox(sim_measurements[end], [9.896156133972513, 0.23998395281014084, 0.7400708222590943]; atol=0.001)
array_isapprox(filtered_beliefs[end].μ, [9.437578760960118, -0.2566452039031286, 0.4989401367606114];
                atol=0.001)
array_isapprox(filtered_beliefs[end].Σ, [0.01682790984107462 5.309144795905572e-5 2.5384310326873065e-5;
                                        5.309144795905572e-5 0.01683608550354201 1.1576863526696874e-7;
                                        2.5384310326873065e-5 1.1576863526696874e-7 0.01683587863740326]; atol=0.001)

@test issymmetric(filtered_beliefs[end].Σ)  
@test isposdef(filtered_beliefs[end].Σ) 

_, info = update_with_info(ekf, b0, action_sequence[1], sim_measurements[1])
@test haskey(info, :innovation_cov)
@test haskey(info, :kalman_gain)
@test haskey(info, :predicted_measurement)
@test size(info.innovation_cov) == (3,3)
@test size(info.kalman_gain) == (3,3) 
@test length(info.predicted_measurement) == 3
