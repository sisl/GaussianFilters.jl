using StaticArrays

@testset "StaticArrays - Linear KF" begin
    rng = StableRNG(0)
    dt = 0.1
    A = @SMatrix [1.0 dt 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 dt; 0.0 0.0 0.0 1.0]
    B = @SMatrix [0.0 0.0; dt 0.0; 0.0 0.0; 0.0 dt]
    W = SMatrix{4,4,Float64}(0.1*I)
    dmodel = LinearDynamicsModel(A, B, W)

    C = @SMatrix [0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0]
    V = SMatrix{2,2,Float64}(0.5*I)
    omodel = LinearObservationModel(C, V)

    kf = KalmanFilter(dmodel, omodel)
    @test kf isa KalmanFilter

    μ0 = @SVector [0.0, 0.0, 0.0, 0.0]
    Σ0 = SMatrix{4,4,Float64}(2.0*I)
    b0 = GaussianBelief(μ0, Σ0)
    @test b0 isa GaussianBelief

    u = @SVector [1.0, 0.5]
    bp = predict(kf, b0, u)
    @test bp isa GaussianBelief

    y = @SVector [0.05, 0.05]
    bn = update(kf, b0, u, y)
    @test bn isa GaussianBelief
end

@testset "StaticArrays - Nonlinear EKF" begin
    rng = StableRNG(0)

    f(x, u) = SVector(x[1] + x[2], x[2] + u[1])
    h(x, u) = SVector(x[1])
    W = SMatrix{2,2,Float64}(0.01*I)
    V = SMatrix{1,1,Float64}(0.1*I)

    dmodel = NonlinearDynamicsModel(f, W)
    omodel = NonlinearObservationModel(h, V)
    ekf = ExtendedKalmanFilter(dmodel, omodel)
    @test ekf isa ExtendedKalmanFilter

    b0 = GaussianBelief(SVector(0.0, 0.0), SMatrix{2,2,Float64}(1.0*I))
    u = SVector(0.1)

    bp = predict(ekf, b0, u)
    @test bp isa GaussianBelief
    @test length(bp.μ) == 2

    y = SVector(0.05)
    bn = update(ekf, b0, u, y)
    @test bn isa GaussianBelief
end

@testset "StaticArrays - UKF" begin
    f(x, u) = SVector(x[1] + u[1], x[2] + u[2])
    h(x, u) = SVector(x[1] + x[2])
    W = SMatrix{2,2,Float64}(0.01*I)
    V = SMatrix{1,1,Float64}(0.1*I)

    dmodel = NonlinearDynamicsModel(f, W)
    omodel = NonlinearObservationModel(h, V)
    ukf = UnscentedKalmanFilter(dmodel, omodel)
    @test ukf isa UnscentedKalmanFilter

    b0 = GaussianBelief(SVector(0.0, 0.0), SMatrix{2,2,Float64}(1.0*I))
    u = SVector(0.1, 0.2)

    bp = predict(ukf, b0, u)
    @test bp isa GaussianBelief

    y = SVector(0.05)
    bn = update(ukf, b0, u, y)
    @test bn isa GaussianBelief
end
