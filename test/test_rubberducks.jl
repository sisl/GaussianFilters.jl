let
    standard_duck_volume = 0.0014 # m^3
    @test isapprox(rubber_ducks_in_earth(0.0014), 3.5769429547897145e21)
end

let
    @test RUBBER_DUCK_PER_SEC == 200
end

let
    @test add_one(1) == 2
end