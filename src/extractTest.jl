include("./classes.jl")
include("./extraction.jl")
include("./plotGMM.jl")

## Main
## 1D example
# Generate GMM
n = 6
w = rand(n)*2
mu = vcat([randn(1)*2 for i=1:n/2],[(randn(1)*3).+10 for i=1:n/2])
sigma = [rand(1,1)*3 for i=1:n]
x = GaussianMixture(w,mu,sigma)

# Extraction
mu_arr = multiple_target_state_extraction(x,0.5)

# Plot
resolution = 1000
x_lim = [-5,15]
plot1dGMM(x, x_lim, resolution)
for mu in mu_arr
	gca()[:axvline](mu,color="k",linestyle="--",alpha=0.5)
end
title("Extraction Test 1D")
xlabel("x")
ylabel("p(x)")
show()


## 2D example
# Generate GMM
n = 6
w = rand(n)
mu = vcat([randn(2) for i=1:n/2],[randn(2).+ [-5,5] for i=1:n/2])
sigma = [Matrix(Diagonal(rand(2))) for i=1:n]
x = GaussianMixture(w, mu, sigma)

# Extraction
mu_arr = multiple_target_state_extraction(x,0.5)

# Plot
xlim = [-10,10]
ylim = [-10,10]
resolution = 100
plot2dGMM(x, xlim, ylim, resolution)
for mu in mu_arr
	plot(mu[1],mu[2],"kx",markeredgewidth=2,markersize=7)
	plot(mu[1],mu[2],"wx",markeredgewidth=1,markersize=6)
end
title("Extraction Test 2D")
xlabel("x_1")
ylabel("x_2")
show()
