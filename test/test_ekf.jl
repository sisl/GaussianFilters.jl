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

# change with algorithm / rng changes

array_isapprox(sim_states[end], [10.7664, 0.213008, -0.463593];atol=0.001)
array_isapprox(sim_measurements[end], [10.9978, 0.491584, -0.33985];atol=0.001)
array_isapprox(filtered_beliefs[end].μ, [10.051, -0.0285158, -0.575077];
                atol=0.001)
array_isapprox(filtered_beliefs[end].Σ, [0.329097 -0.00281178 0.000218104;
                                        -0.00281178 0.0168616 -1.87166e-6;
                                        0.000218104 -1.87166e-6 0.0168371]; atol=0.001)
