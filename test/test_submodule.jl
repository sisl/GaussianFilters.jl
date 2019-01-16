let
    # Test first function here
    r1, r2 = your_function(3, 1.5)
    @test r1 == 4.5
    @test r2 == 1.5

    # Test corner case to function here
    r1, r2 = your_function(1, -42)
    @test r1 == 0
    @test r2 == 0
end

let
    # Here's an example of using the function array_isapprox 
    array = [1.0 1.001 0.999]
    array_isapprox(array, 1.0, atol=0.002)
end