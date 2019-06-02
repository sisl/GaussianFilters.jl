using PyPlot
using LinearAlgebra
using Distributions

"""
Plot a 2D GMM

Arguments:
- `x::GaussianMixture` GMM State distribution. [Gaussian Mixture]
- `xlim::Vector` [xmin,xmax] for plotting
- `ylim::Vector` [ymin,ymax] for plotting
- `Resolution::Integer` Number of points between xlim (and ylim)
"""
function plot2dGMM(x::GaussianMixture, xlim, ylim, resolution)

	# Form grid
	x_grid = range(xlim[1],stop=xlim[2],length=resolution)
	y_grid = range(ylim[1],stop=ylim[2],length=resolution)
	grid = [[x_grid[i],y_grid[j]] for i = 1:length(x_grid) for j = 1:length(y_grid)]
	grid = hcat(grid...)
	surf = zeros(resolution,resolution)

	# Build surface
	N = x.N
	for i = 1:N
	    w = x.w[i]
	    μ = x.μ[i]
	    Σ = x.Σ[i]

	    surf_i = w.*pdf(MvNormal(μ,Σ), grid)
	    surf_i = reshape(surf_i,resolution,resolution)
	    surf = surf + surf_i
	end

	figure()
	imshow(surf,origin="lower",extent=vcat(xlim,ylim))
	xlabel("x_1")
	ylabel("x_2")
	colorbar()
end


"""
Plot a 1D GMM

Arguments:
- `x::GaussianMixture` GMM State distribution. [Gaussian Mixture]
- `xlim::Vector` [xmin,xmax] for plotting
- `Resolution::Integer` Number of points between xlim
"""
function plot1dGMM(x::GaussianMixture, xlim, resolution)

	# Form grid
	x_grid = range(xlim[1],stop=xlim[2],length=resolution)
	y = zeros(resolution)

	# Build surface
	N = x.N
	for i = 1:N
	    w = x.w[i]
	    μ = x.μ[i][1]
	    Σ = x.Σ[i][1]

	    sigma = sqrt(Σ)
	    y_i = w.*pdf.(Normal(μ,sigma), x_grid)
	    y = y + y_i
	end

	figure()
	plot(x_grid,y)
	xlabel("x")
	ylabel("p(x)")
end
