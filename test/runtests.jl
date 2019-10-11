# Packages required for testing
using Test
using Random
using Logging
using LinearAlgebra
using Distributions
using NBInclude

# Package Under Test
using GaussianFilters

# Set logging level
global_logger(SimpleLogger(stderr, Logging.Debug))

# Fix randomness during tests
Random.seed!(0)

# Check equality of two arrays
@inline function array_isapprox(x::AbstractArray{F},
                  y::AbstractArray{F};
                  rtol::F=sqrt(eps(F)),
                  atol::F=zero(F)) where {F<:AbstractFloat}

    # Easy check on matching size
    if length(x) != length(y)
        return false
    end

    for (a,b) in zip(x,y)
        @test isapprox(a,b, rtol=rtol, atol=atol)
    end
end

# Check if array equals a single value
@inline function array_isapprox(x::AbstractArray{F},
                  y::F;
                  rtol::F=sqrt(eps(F)),
                  atol::F=zero(F)) where {F<:AbstractFloat}

    for a in x
        @test isapprox(a, y, rtol=rtol, atol=atol)
    end
end

# Define package tests

@time @testset "GaussianFilter Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")

    @time @testset "Models Tests" begin
        include(joinpath(testdir, "test_models.jl"))
    end

    @time @testset "GaussianFilter GM-PHD Testing" begin
        include(joinpath(testdir, "test_gmphd.jl"))
    end

    @time @testset "GaussianFilter KF Testing" begin
        include(joinpath(testdir, "test_kf.jl"))
    end

    @time @testset "GaussianFilter EKF Testing" begin
        include(joinpath(testdir, "test_ekf.jl"))
    end

    @time @testset "GaussianFilter UKF Testing" begin
        include(joinpath(testdir, "test_ukf.jl"))
    end

    @testset "Notebooks testing" begin 
        nbdir = joinpath(dirname(@__DIR__), "notebooks")
        for d in readdir(nbdir)
            if endswith(d, ".ipynb")
                path = joinpath(nbdir, d)
                @testset "$d" begin
                    @nbinclude path
                end
            end
        end
    end

end
