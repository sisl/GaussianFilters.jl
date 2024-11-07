Random.seed!(2019)

# Dynamics
dt=0.1
function dynamics(x,u)
    v,ϕ = u
    xp = x
    xp += dt*([v*cos(x[3]);v*sin(x[3]); ϕ])
    xp[3] = xp[3]
    return xp
end

W = dt*0.05*Matrix{Float64}(I,3,3)

dmodel = NonlinearDynamicsModel(dynamics,W)

# Measurement
observe(x,u) = [norm(x[1:2])]

V = 0.1*Matrix{Float64}(I,1,1)

omodel = NonlinearObservationModel(observe,V)

# Build UKF
ukf = UnscentedKalmanFilter(dmodel,omodel);

# Simulation
times = 0:dt:5
action_sequence = [[1.0, 2*sin(pi/2*t)] for t in times]

b0 = GaussianBelief([0.0,0.0,0.0], 0.01*Matrix{Float64}(I,3,3))

sim_states, sim_measurements = simulation(ukf,b0,action_sequence);

# Filter
filtered_beliefs = run_filter(ukf, b0, action_sequence, sim_measurements)
μ, Σ = unpack(filtered_beliefs)

# Tests
@test ukf.d isa NonlinearDynamicsModel
@test ukf.o isa NonlinearObservationModel
@test length(sim_measurements) == length(action_sequence)

points, w_μ, w_Σ = @inferred unscented_transform(b0, ukf.λ, ukf.α, ukf.β)
bp = @inferred unscented_transform_inverse(points, w_μ, w_Σ)

first_update = update(ukf, b0, action_sequence[1], sim_measurements[1])
first_predict = predict(ukf, b0, action_sequence[1])
first_predict_measure = measure(ukf, first_predict, sim_measurements[1],
                        u = action_sequence[1])
array_isapprox(first_update.μ, first_predict_measure.μ)
array_isapprox(first_update.Σ, first_predict_measure.Σ)


@test length(filtered_beliefs[end].μ) == length(b0.μ)
@test size(filtered_beliefs[end].Σ) == size(b0.Σ)
@test size(μ) == (length(filtered_beliefs), length(filtered_beliefs[1].μ))
@test size(Σ) == (length(filtered_beliefs), size(filtered_beliefs[1].Σ)...)

# Test mathematical properties
@test all(isfinite.(sim_states[end]))  
@test all(isfinite.(sim_measurements[end])) 
@test length(sim_measurements[end]) == 1  
@test length(sim_states[end]) == 3 

# Test filtered belief properties
@test all(isfinite.(filtered_beliefs[end].μ))  
@test all(isfinite.(filtered_beliefs[end].Σ)) 
@test issymmetric(filtered_beliefs[end].Σ)  
@test isposdef(filtered_beliefs[end].Σ) 
@test tr(filtered_beliefs[end].Σ) < 2.0 
