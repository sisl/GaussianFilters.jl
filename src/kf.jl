"""
    Work into Kalman update, separate by predict, measure,
    update where update = predict, measure
"""
function KalmanUpdate(μ, Σ, u, y, A, B, C, D, W, V)
    μ_pred = A*μ + B*u
    Σ_pred = A*Σ*A' + W

    K = Σ_pred*C'*inv(C*Σ_pred*C'+V)
    μ = μ_pred + K*(y-(C*μ_pred + D*u))
    Σ = (I - K*C)*Σ_pred
    return μ, Σ
end
