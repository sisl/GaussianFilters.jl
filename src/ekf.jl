"""
    Work into extended Kalman update, separate by predict, measure,
    update where update = predict, measure, use autodiff to compute G and H
"""

function ExtendedKalmanUpdate(μ, Σ, u, y, g, h, W, V;dt=0.001)
    μ_pred = g(μ, u;dt=dt)
    # G = dg/dmu
    Σ_pred = G*Σ*G' + W

    # H = dh/dmu
    K = Σ_pred*H'*inv(H*Σ_pred*H'+V)
    μ = μ_pred + K*(y-h(μ_pred))
    Σ = (I - K*H)*Σ_pred
    return μ, Σ
end
