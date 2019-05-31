include("./classes.jl")
include("./prune.jl")
include("./plotGMM.jl")
using Printf

## Main
# 1d example
n = 1000

w = rand(n)
mu = vcat([randn(1)*3 for i=1:n/2],[(randn(1)*3).+10 for i=1:n/2])
sigma = [rand(1,1) for i=1:n]

x = GaussianMixture(w,mu,sigma)

T = 0.0			# threshold
U = 10.0		# clustering threshold
J_max = 100

x_new = prune(x,T,U,J_max)

xlim = [-5,15]
resolution = 1000
plot1dGMM(x, xlim, resolution)
title_str = @sprintf "Before Pruning: N=%d" n
title(title_str)

plot1dGMM(x_new, xlim, resolution)
title_str = @sprintf("After Pruning: N=%d", length(x_new.w))
title(title_str)
show()
