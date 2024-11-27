# Arbitrary GM-PHD example
let
    rng = StableRNG(0)

    x0 = GaussianMixture([1.0],[randn(rng, 2)],[Matrix(1.0I,2,2)])
    Z = [randn(rng, 1)]
    k(z) = 0
    gamma = GaussianMixture([1.0],[randn(rng, 2)],[Matrix(1.0I,2,2)])
    sp_dyn = [Dynamics(Matrix(1.0I,2,2),Matrix(1.0I,2,2),[0.0, 0.0])]
    sp = Spawn(GaussianMixture([1.0],[randn(rng, 2)],[Matrix(1.0I,2,2)]), sp_dyn)
    dyn = Dynamics(Matrix(1.0I,2,2),Matrix(1.0I,2,2))
    mea = Measurement(randn(rng, 1,2),ones(1,1))

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
    array_isapprox(s.w, [0.5, 0.5, 0.25, 0.19792404040719125, 0.5347173063952059, 0.26735865319760294], atol=0.002)
    @test length(s.Σ) == s.N
    @test length(s.μ) == s.N

    # Prune and test pruning
    T = 10^-5 # Truncation threshold
    U = 0.1 # Merging threshold
    J_max = 10; # Max number of Gaussian terms
    p = prune(s, T, U, J_max)
    @test p.N == 4
    array_isapprox(p.w, [0.8020759595928089, 0.5, 0.75, 0.19792404040719125], atol=0.002)
    @test length(p.μ) == p.N

    # Test single step updating
    sp = update(phd1, x0, Z, T, U, J_max)
    @test sp.w == p.w
    @test sp.μ == p.μ

    # Test state extraction at different thresholds
    states1 = multiple_target_state_extraction(sp, 0.5)
    @test length(states1) == 2
    array_isapprox(states1[1], [-0.16045023924680923, 0.4344660539101717], atol=0.002)
    states2 = multiple_target_state_extraction(sp, 2.0)
    @test length(states2) == 0
    array_isapprox(states1[2], [-0.3811118617867437, -0.020452508120591215], atol=0.002)
end

# Surveillance example
let
    rng = StableRNG(0)
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
        x_new = [  Dyns.A*xsim[t][i] + rand(rng, MVD,1)[:] for i=1:length(xsim[t]) ]
        if t == 13
            MVDspawn = MvNormal(Q_β)
            x1 = xsim[t][1]
            x_spawn = Dyns.A*x1 + rand(rng, MVDspawn,1)[:]
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
        z = [ Meas.C*xsim[t][i] + rand(rng, MVD,1)[:] for i=1:length(xsim[t])]
        x_new_pruned = update(phd,x[t],z,T,U,J_max)
        push!(x, x_new_pruned)
    end

    sf = x[end]
    @test sf.N == 6
    array_isapprox(sf.w, [0.9859331369043042, 0.974981461116766, 0.9417729637426844,
                          0.9135806098534265, 0.020487604998975638, 0.020487604998975638], atol=0.002)
    @test length(sf.Σ) == 6

    esf = multiple_target_state_extraction(sf,0.5)
    @test length(esf) == 4
    array_isapprox(esf[1], [-252.02649932803052, -359.93362714587613, 4.6283940976004185, 8.39795059665333], atol=0.002)
    array_isapprox(esf[2], [-259.95843400206553, -222.24691047655602, -3.2305683889212893, 1.5813642133068648], atol=0.002)
    array_isapprox(esf[3], [37.31131019093644, 308.138415366749, -25.376170638745595, -0.49215530411660485], atol=0.002)
    array_isapprox(esf[4], [-287.2383000244064, -489.18789830748494, -14.57194769862216, -10.073401647546325], atol=0.002)
end
