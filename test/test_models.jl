@testset "linear" begin
    # test linear 
    # x dim = 4, u dim = 2
    A = 1.0 * Matrix(I, (4, 4))
    B = [1. 0.; 0. 1.;0. 0.; 0. 0.]
    C = 1.0 * Matrix(I, (3, 4))
    D = zeros(3, 2)
    W = 0.5 * Matrix(I, (4, 4))
    V = 0.1 * Matrix(I, (3, 3))
    dm = LinearDynamicsModel(A, B, W)
    om = LinearObservationModel(C, D, V)

    rng = MersenneTwister(1)

    x0 = ones(4)
    u0 = ones(2)

    x1 = predict(dm, x0, u0)
    @test x1 == [2.0, 2.0, 1.0, 1.0]
    @inferred predict(dm, x0, u0)

    x2 = predict(dm, x0, u0, rng)
    @test x2[1] > 1.0

    o1 = measure(om, x0, u0)
    @inferred measure(om, x0, u0)
    @test o1 == [1.0, 1.0, 1.0]
    o2 = measure(om, x0, u0, rng)
end

# test non linear
@testset "non linear" begin 

    rng = MersenneTwister(1)

    W = 0.001*Matrix{Float64}(I,4,4)
    V = 0.3*Matrix{Float64}(I,4,4)

    dm = NonlinearDynamicsModel((x, u) -> x.^2 + u, W);
    om = NonlinearObservationModel((x, u) -> x.^4, V)


    x0 = ones(4)
    u0 = ones(4)

    x1 = predict(dm, x0, u0)
    @inferred predict(dm, x0, u0)
    @test x1 == 2*ones(4)
    x2 = predict(dm, x1, u0, rng) 
    @test all(i -> i != 5.0, x2)

    o1 = measure(om, x0, u0)
    @inferred measure(om, x0, u0)
    @test o1 == x0
    o2 = measure(om, x0, u0, rng)
    @test all(x -> !isapprox(x, 1.0), o2)
end
