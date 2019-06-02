include("./classes.jl")
include("./extraction.jl")
include("./plotGMM.jl")

## Main
# 1D example
n = 6

w = rand(n)*2
mu = vcat([randn(1)*2 for i=1:n/2],[(randn(1)*3).+10 for i=1:n/2])
sigma = [rand(1,1)*3 for i=1:n]

x = GaussianMixture(w,mu,sigma)
mu_arr = multiple_target_state_extraction(x,0.5)

T = 0.0			# threshold
U = 10.0		# clustering threshold
J_max = 100

resolution = 1000
x_lim = [-5,15]
plot1dGMM(x, x_lim, resolution)
for mu in mu_arr
	gca()[:axvline](mu,color="k",linestyle="--",alpha=0.5)
end
title("Extraction Test 1D")
show()


# 2D example
n = 6

w = rand(n)
mu = vcat([randn(2) for i=1:n/2],[randn(2).+ [-5,5] for i=1:n/2])
sigma = [Matrix(Diagonal(rand(2))) for i=1:n]

T = 0.0			# threshold
U = 10.0		# clustering threshold
J_max = 100		# maximum number of points

x = GaussianMixture(w, mu, sigma)
mu_arr = multiple_target_state_extraction(x,0.5)

xlim = [-10,10]
ylim = [-10,10]
resolution = 100
plot2dGMM(x, xlim, ylim, resolution)
for mu in mu_arr
	plot(mu[1],mu[2],"kx",markeredgewidth=2,markersize=7)
	plot(mu[1],mu[2],"wx",markeredgewidth=1,markersize=6)
end
title("Extraction Test 2D")

show()
