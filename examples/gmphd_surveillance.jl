using GaussianFilters

using LinearAlgebra
using Distributions

Δ = 1 # Sampling Period
σ_p = 5


F = [ Matrix{Float64}(I,2,2) Δ*Matrix{Float64}(I,2,2) ;
    zeros(Float64,2,2) Matrix{Float64}(I,2,2) ]

Q  = σ_p^2 * [ (Δ^4)/4*Matrix{Float64}(I,2,2) (Δ^3)/4*Matrix{Float64}(I,2,2) ;
            (Δ^3)/4*Matrix{Float64}(I,2,2) (Δ^2)*Matrix{Float64}(I,2,2) ]

Dyns = Dynamics(F,Q);

σ_m = 10

C = [ Matrix{Float64}(I,2,2) zeros(Float64,2,2) ]
R = σ_m^2 * Matrix{Float64}(I,2,2)

Meas = Measurement(C,R);

w1 = 0.1
w2 = 0.1
m_γ1 = [ 250 , 250 , 0 , 0]
m_γ2 = [ -250 , -250 , 0 , 0 ]
P_γ = Matrix(Diagonal([100 , 100 , 25 , 25]));

Q_β = Matrix(Diagonal([100.0 , 100.0 , 400.0 , 400.0]))

β = GaussianMixture([w1, w2] , [m_γ1, m_γ2] , [P_γ, P_γ])
spawn = Spawn(β , [Dyns,Dyns]);

# Creates a uniform distriubtion given a single measurement z [px,py]
function κ(z)
    V = 4*10^6 #Surveillance volume
    λ_c = 12.5*10^-6 #Average clutter returns per unit volume
    return λ_c*V*(1/(2000^2))
end;

Ps = 0.99 # Probability of survival
Pd = 0.98 # Probability of detection

T = 10^-5 # Truncation threshold
U = 4 # Merging threshold
J_max = 10; # Max number of Gaussian terms

w_target1 = 1;
w_target2 = 1;
mu_target1 = [ -300.0 , -300.0 , 0.0 , 0.0 ]
mu_target2 = [ 300.0 , 300.0 , 0.0 , 0.0 ]
P_agent1 = Matrix(Diagonal([100.0 , 100.0 , 25.0 , 25.0]))
P_agent2 = P_agent1
γ = GaussianMixture([w_target1, w_target2] , [mu_target1, mu_target2] ,
                    [P_agent1, P_agent2]);

xsim = [γ.μ]
zs = []

for t = 1 : Δ : 100
    # Simulate state
    MVD = MvNormal(Dyns.Q)
    x_new = [  Dyns.A*xsim[t][i] + rand(MVD,1)[:] for i=1:length(xsim[t]) ]
    if t == 66
        MVDspawn = MvNormal(Q_β)
        x1 = xsim[t][1]
        x_spawn = Dyns.A*x1 + rand(MVDspawn,1)[:]
        push!(x_new,x_spawn)
    end
    if t == 40
        x_birth = m_γ2
        push!(x_new,x_birth)
    end
    push!(xsim,x_new)
    
    # Simulate Measurements
    MVD = MvNormal(Meas.R)
    z = [ Meas.C*xsim[t][i] + rand(MVD,1)[:] for i=1:length(xsim[t])]
    push!(zs, z)
end

# Instantiate the GM-PHD filter
phd = PHDFilter(γ,spawn,Dyns,Meas,Ps,Pd,κ)

# Initialize state of objects: x = γ
x =  [γ]

# Run PHD filter
for t = 1:Δ:100
    x_new_pruned = update(phd, x[t], zs[t] ,T,U,J_max)
    push!(x, x_new_pruned)
end

using Plots

p1 = plot(legend=false)
p2 = plot(legend=false)

for (t,x) in enumerate(x)
    mu_arr = multiple_target_state_extraction(x,0.5)
    for mu in mu_arr
        scatter!(p1, [t], [mu[1]],color=:black)
        scatter!(p2, [t],[mu[2]],color=:black)
    end
    scatter!(p1, [t],[xsim[t][1][1]],color=:blue)
    scatter!(p1, [t],[xsim[t][2][1]],color=:blue)
    if t > 66
        scatter!(p1, [t],[xsim[t][4][1]],color=:blue)
    end
    if t > 40
        scatter!(p1, [t],[xsim[t][3][1]],color=:red)
    end
    scatter!(p2, [t],[xsim[t][1][2]],color=:blue)
    scatter!(p2, [t],[xsim[t][2][2]],color=:blue)
    if t > 66
        scatter!(p2, [t], [xsim[t][4][2]],color=:blue)
    end
    if t > 40
        scatter!(p2, [t], [xsim[t][3][2]], color=:blue)
    end
end
xlabel!(p1, "time step")
ylabel!(p1, "x-coordinate")
xlabel!(p2, "time step")
ylabel!(p2, "y-coordinate")
plot(p1, p2, layout=(2,1))
