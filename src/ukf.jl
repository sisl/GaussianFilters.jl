"""
Unscented Transform

Convert a single mean and covariance to set of 2n+1 sigma points
with n being the dimensionality of the state space.
Uses formulation from ProbRob for α/β parameters on a separate covariance
weighting, although this is not necessary.

"""
function UT(μ, Σ, λ=2, α=1, β=0)
    n = length(μ)

    # compute weights
    w_m = 1/(2.0*(n+λ))*ones(2*n+1,1)
    w_m[1] = λ/(λ+n)
    w_c = copy(w_m)
    w_c[1] += (1 - α^2 + β) # Per ProbRob formulation

    s = cholesky(Σ).L
    points = zeros(n,2*n+1)
    points[:,1] = μ
    for i in 1:n
        points[:,2*i] = μ + sqrt(n+λ)*s[:,i]
        points[:,2*i+1] = μ - sqrt(n+λ)*s[:,i]
    end

    return points, w_m, w_c
end

"""
Inverse Unscented Transform

Convert a 2n+1 sigma points and weights back to a single measure for
mean and covariance.

Uses formulation from ProbRob for α/β parameters on a separate covariance
weighting, although this is not necessary.

"""
function invUT(points, w_m, w_c)
    μ = points*w_m
    diff = points .- μ
    scaled = diff.*sqrt.(w_c)'
    Σ = scaled*scaled'

    return μ, Σ
end

"""
Unscented Kalman Update

Rework into separate predict and measure functions, with update=predict+measure
"""
function UnscentedKalmanUpdate(μ, Σ, u, y, g, h, W, V ;dt=0.001)
    n = length(μ)

    # predict
    points, w_m, w_c = UT(μ,Σ)
    pointsp = zeros(size(points))
    for i in 1:(2*n+1)
        pointsp[:,i] = g(points[:,i], u;dt=dt)
    end
    μ_pred, Σb_pred = invUT(pointsp, w_m, w_c)
    Σ_pred = Σb_pred + W

    # update
    points, w_m, w_c = UT(μ_pred,Σ_pred)
    y_preds = zeros(length(y),(2*n+1))
    for i in 1:(2*n+1)
        y_preds[:,i] = h(points[:,i])
    end
    y_pred = y_preds*w_m

    scaled_ydiff = (y_preds .- y_pred).*sqrt.(w_c)'
    Σ_Y = scaled_ydiff*scaled_ydiff' + V

    scaled_xdiff = (points .- μ_pred).*sqrt.(w_c)'
    Σ_XY = scaled_xdiff*scaled_ydiff'

    μ = μ_pred + Σ_XY*inv(Σ_Y)*(y-y_pred)
    Σ = Σ_pred - Σ_XY*inv(Σ_Y)*Σ_XY'
    Σ = Symmetric(Σ)
    return μ, Σ
end
