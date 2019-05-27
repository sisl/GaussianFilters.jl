include("./classes.jl")
include("./prune.jl")
using PyPlot


## Main
# function GaussianMixture(w, μ::Vector{Vector{T}}, Σ) where T
n = 1000
x = GaussianMixture(rand(n),vcat([randn(1) for i=1:n/2],[randn(1).+10 for i=1:n/2]),[rand(1,1) for i=1:n] )

T = 0.1		# threshold
U = 100.0		# clustering threshold
J_max = 5

x_new = prune(x,T,U,J_max)

figure()
plot(x.μ,x.w,".")
xlabel("μ")
ylabel("w")

figure()
plot(x_new.μ,x_new.w,".")
xlabel("μ")
ylabel("w")

show()
