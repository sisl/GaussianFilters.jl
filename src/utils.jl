"""
    belief_ellipse(b::GaussianBelief, P::Float=0.95; δ::Number=5)

Construct and return the x and y points of a 2D gaussian belief,
with P being the total probability captured by the ellipse (P ∈ (0,1)),
and δ the degree increment between points.
"""
function belief_ellipse(b::GaussianBelief, P::Number=0.95; δ::Number=5)
    @assert 0<P<1
    @assert 0<δ<360

    θ = (pi/180)*collect(0:δ:360)
    rad_w = (-2*log(1-P))^0.5

    w1 = rad_w*cos.(θ)
    w2 = rad_w*sin.(θ)

    Sig12 = cholesky(b.Σ)
    x=Sig12.L*[w1 w2]' .+ b.μ
    return x[1,:], x[2,:]
end

"""
    belief_ellipse(μ::AbstractVector, Σ::AbstractMatrix, P::Float=0.95; δ::Number=5)

Construct and return the x and y points of a 2D gaussian belief,
with P being the total probability captured by the ellipse (P ∈ (0,1)),
and δ the degree increment between points.
"""

function belief_ellipse(μ::AbstractVector, Σ::AbstractMatrix, P::Number=0.95; δ::Number=5)
    @assert 0<P<1
    @assert 0<δ<360

    θ = (pi/180)*collect(0:δ:360)
    rad_w = (-2*log(1-P))^0.5

    w1 = rad_w*cos.(θ)
    w2 = rad_w*sin.(θ)

    Sig12 = cholesky(Σ)
    x = Sig12.L*[w1 w2]' .+ μ
    return x[1,:], x[2,:]
end
