using POMDPs
using Distributions: MvNormal

@testset "POMDPs.update extension" begin
    dt = 0.1
    A = [1.0 dt; 0.0 1.0]
    B = [0.0; dt][:, :]
    W = 0.01 * Matrix{Float64}(I, 2, 2)
    dmodel = LinearDynamicsModel(A, B, W)

    C = [1.0 0.0]
    V = 0.1 * Matrix{Float64}(I, 1, 1)
    omodel = LinearObservationModel(C, V)

    kf = KalmanFilter(dmodel, omodel)
    b0 = GaussianBelief([0.0, 0.0], Matrix{Float64}(I, 2, 2))

    # POMDPs.update should produce the same result as GaussianFilters.update
    a = [1.0]
    o = [0.5]
    b_pomdps = POMDPs.update(kf, b0, a, o)
    b_native = GaussianFilters.update(kf, b0, a, o)
    @test b_pomdps.μ == b_native.μ
    @test b_pomdps.Σ == b_native.Σ

    # initialize_belief should be the identity on a GaussianBelief
    @test POMDPs.initialize_belief(kf, b0) === b0

    # initialize_belief from an MvNormal extracts mean and covariance
    d = MvNormal([1.0, 2.0], [4.0 0.0; 0.0 9.0])
    b_from_d = POMDPs.initialize_belief(kf, d)
    @test b_from_d isa GaussianBelief
    @test b_from_d.μ == [1.0, 2.0]
    @test b_from_d.Σ == [4.0 0.0; 0.0 9.0]
end
