using LinearAlgebra
using Distributions

include("./classes.jl")
include("./PHDStep.jl")

σ_p = 5 #Standard deviation - process noise
σ_m = 10 #Standard deviation - measurement noise
Δ = 1 #Sampling Period

# Dynamics Model - Gaussian
F = [ Matrix{Float64}(I,2,2) Δ*Matrix{Float64}(I,2,2) ;
    zeros(Float64,2,2) Matrix{Float64}(I,2,2) ]

Q  = σ_p^2 * [ (Δ^4)/4*Matrix{Float64}(I,2,2) (Δ^3)/2*Matrix{Float64}(I,2,2) ;
            (Δ^3)/2*Matrix{Float64}(I,2,2) (Δ^2)*Matrix{Float64}(I,2,2) ]

Dyns = Dynamics(F,Q)

# Measurement Model
C = [ Matrix{Float64}(I,2,2) zeros(Float64,2,2) ]
R = σ_m^2 * Matrix{Float64}(I,2,2)

Meas = Measurement(C,R)

# Spontaneous birth model
w1 = 0.1
w2 = 0.1
m_γ1 = [ 250 , 250 , 0 , 0]
m_γ2 = [ -250 , -250 , 0 , 0 ]
P_γ = Matrix(Diagonal([100 , 100 , 25 , 25]))

# Spawn (from existing targets) model
Q_β = Matrix(Diagonal([100 , 100 , 400 , 400]))

β = GaussianMixture([w1, w2] , [m_γ1, m_γ2] , [P_γ, P_γ])
spawn = Spawn(β , [Dyns,Dyns])

# Clutter model
""" Creates a uniform distriubtion given a single measurement z [px,py] """
function κ(z)
    V = 4*10^6 #Surveillance volume
    λ_c = 12.5*10^-6 #Average clutter returns per unit volume
    #κ
    λ_c*V*(1/(2000^2))
end

# PHD parameters
T = 10^-5 #Truncation threshold
U = 4 #Merging threshold
J_max = 100 #Max number of Gaussian terms
Ps = 0.99 #Probability of survival
Pd = 0.98 #Probability of detection

# Start with two existing targets
w_target1 = 0.2;
w_target2 = 0.2;
mu_target1 = [ -300.0 , -300.0 , 0.0 , 0.0 ]
mu_target2 = [ 300.0 , 300.0 , 0.0 , 0.0 ]
P_agent1 = Matrix(Diagonal([100.0 , 100.0 , 25.0 , 25.0]))
P_agent2 = P_agent1
γ = GaussianMixture([w_target1, w_target2] , [mu_target1, mu_target2] ,
                    [P_agent1, P_agent2])

# Build PHD
phd = PHDFilter(γ,spawn,Dyns,Meas,Ps,Pd,κ)


# Simulate system with PHD filter
# Initilaise s.t. filter knows state of objects: x = γ
w = []
mu = [[]]
P = [zeros(1,1)]

#x = Array{GaussianMixture}
x =  [γ]
for t = 1:Δ:100
    MVD = MvNormal(Meas.R)
    z = [ Meas.C*x[t].μ[i] + rand(MVD,1)[:] for i=1:x[t].N]
    x_new = step(x[t],z,phd)
    println(typeof(x[1]))
    println(typeof(x_new))
    push!(x, x_new)
end
