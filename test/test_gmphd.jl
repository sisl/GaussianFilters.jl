# Arbitrary GM-PHD example
let
    Random.seed!(0) # reset seed

    x0 = GaussianMixture([1.0],[randn(2)],[Matrix(1.0I,2,2)])
    Z = [randn(1)]
    k(z) = 0
    gamma = GaussianMixture([1.0],[randn(2)],[Matrix(1.0I,2,2)])
    sp_dyn = [Dynamics(Matrix(1.0I,2,2),Matrix(1.0I,2,2),[0.0, 0.0])]
    sp = Spawn(GaussianMixture([1.0],[randn(2)],[Matrix(1.0I,2,2)]), sp_dyn)
    dyn = Dynamics(Matrix(1.0I,2,2),Matrix(1.0I,2,2))
    mea = Measurement(randn(1,2),ones(1,1))

    # Should construct identically
    phd1 = PHDFilter(gamma, sp, dyn, mea, 0.5, 0.5, k)
    phd2 = PHDFilter(gamma, sp, [dyn], mea, 0.5, 0.5, k) # should be same
    @test phd1.γ.N == 1
    @test phd2.γ.N == 1
    @test phd1.dyn == phd2.dyn

    # Predict/Measure and test update
    xp = predict(phd1, x0)
    s = measure(phd1, xp, Z)
    @test s.N == 6
    array_isapprox(s.w, [0.5, 0.5, 0.25, 0.425012, 0.383326, 0.191663], atol=0.002)
    @test length(s.Σ) == s.N
    @test length(s.μ) == s.N

    # Prune and test pruning
    T = 10^-5 # Truncation threshold
    U = 0.1 # Merging threshold
    J_max = 10; # Max number of Gaussian terms
    p = prune(s, T, U, J_max)
    @test p.N == 2
    array_isapprox(p.w, [0.925012, 1.32499], atol=0.002)
    @test length(p.μ) == p.N

    # Test single step updating
    sp = update(phd1, x0, Z, T, U, J_max)
    @test sp.w == p.w
    @test sp.μ == p.μ

    # Test state extraction at different thresholds
    states1 = multiple_target_state_extraction(sp, 0.8)
    @test length(states1) == 2
    array_isapprox(states1[1], [-0.132261, 0.598847], atol=0.002)

    states2 = multiple_target_state_extraction(sp, 2.0)
    @test length(states2) == 0
end

# Surveillance example
let
    Random.seed!(0)

    # Dynamics
    Δ = 1 # Sampling Period
    σ_p = 5
    F = [ Matrix{Float64}(I,2,2) Δ*Matrix{Float64}(I,2,2) ;
        zeros(Float64,2,2) Matrix{Float64}(I,2,2) ]
    Q  = σ_p^2 * [ (Δ^4)/4*Matrix{Float64}(I,2,2) (Δ^3)/4*Matrix{Float64}(I,2,2) ;
                (Δ^3)/4*Matrix{Float64}(I,2,2) (Δ^2)*Matrix{Float64}(I,2,2) ]
    Dyns = Dynamics(F,Q);

    # Measurement
    σ_m = 10
    C = [ Matrix{Float64}(I,2,2) zeros(Float64,2,2) ]
    R = σ_m^2 * Matrix{Float64}(I,2,2)
    Meas = Measurement(C,R);

    # Birth
    w1 = 0.1
    w2 = 0.1
    m_γ1 = [ 250 , 250 , 0 , 0]
    m_γ2 = [ -250 , -250 , 0 , 0 ]
    P_γ = Matrix(Diagonal([100 , 100 , 25 , 25]));

    # Spawn
    Q_β = Matrix(Diagonal([100.0 , 100.0 , 400.0 , 400.0]))
    β = GaussianMixture([w1, w2] , [m_γ1, m_γ2] , [P_γ, P_γ])
    spawn = Spawn(β , [Dyns,Dyns]);

    # Clutter
    function κ(z)
        V = 4*10^6 #Surveillance volume
        λ_c = 12.5*10^-6 #Average clutter returns per unit volume
        return λ_c*V*(1/(2000^2))
    end;

    # Misc
    Ps = 0.99 # Probability of survival
    Pd = 0.98 # Probability of detection
    T = 10^-5 # Truncation threshold
    U = 4 # Merging threshold
    J_max = 10; # Max number of Gaussian terms

    # Instantiation
    w_target1 = 1;
    w_target2 = 1;
    mu_target1 = [ -300.0 , -300.0 , 0.0 , 0.0 ]
    mu_target2 = [ 300.0 , 300.0 , 0.0 , 0.0 ]
    P_agent1 = Matrix(Diagonal([100.0 , 100.0 , 25.0 , 25.0]))
    P_agent2 = P_agent1
    γ = GaussianMixture([w_target1, w_target2] , [mu_target1, mu_target2] ,
                        [P_agent1, P_agent2]);

    # Simulation
    xsim = [γ.μ]
    Tf = 20
    for t = 1 : Δ : Tf
        MVD = MvNormal(Dyns.Q)
        x_new = [  Dyns.A*xsim[t][i] + rand(MVD,1)[:] for i=1:length(xsim[t]) ]
        if t == 13
            MVDspawn = MvNormal(Q_β)
            x1 = xsim[t][1]
            x_spawn = Dyns.A*x1 + rand(MVDspawn,1)[:]
            push!(x_new,x_spawn)
        end
        if t == 8
            x_birth = m_γ2
            push!(x_new,x_birth)
        end
        push!(xsim,x_new)
    end

    # Filter
    phd = PHDFilter(γ,spawn,Dyns,Meas,Ps,Pd,κ)
    x =  [γ]
    for t = 1:Δ:Tf
        MVD = MvNormal(Meas.R)
        z = [ Meas.C*xsim[t][i] + rand(MVD,1)[:] for i=1:length(xsim[t])]
        x_new_pruned = update(phd,x[t],z,T,U,J_max)
        push!(x, x_new_pruned)
    end

    sf = x[end]
    @test sf.N == 6
    array_isapprox(sf.w, [0.999298, 0.983239, 0.95002,
                        1.03055, 0.0204876, 0.000549536], atol=0.002)
    @test length(sf.Σ) == 6

    esf = multiple_target_state_extraction(sf,0.5)
    @test length(esf) == 4
    array_isapprox(esf[1], [-282.13, -724.875, -1.9644, -46.7133], atol=0.002)
    array_isapprox(esf[4], [-313.933, -301.799, -3.09298, -1.9735], atol=0.002)
end
