# Packages required for testing
using Test
using Random
using Logging

# Package Under Test
using JuliaPackageTemplate

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
@time @testset "JuliaPackageTemplate Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    @time @testset "JuliaPackageTemplate.YourSubmodule" begin
        include(joinpath(testdir, "test_submodule.jl"))
    end
    @time @testset "JuliaPackageTemplate.SingleSatellite" begin
        include(joinpath(testdir, "test_rubberducks.jl"))
    end
end