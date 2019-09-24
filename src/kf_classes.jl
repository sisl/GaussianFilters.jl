### Motion/Measurement Models ###

"""
`DynamicsModel` is an abstract type to encapsulate linear
and nonlinear dynamics models
"""
abstract type DynamicsModel end

"""
`ObservationModel` is an abstract type to encapsulate linear
and nonlinear observation models
"""
abstract type ObservationModel end

"""
    LinearDynamicsModel(A::AbstractMatrix,B::AbstractMatrix,W::Symmetric)
    LinearDynamicsModel(A::AbstractMatrix,B::AbstractMatrix,W::AbstractMatrix)

Construct linear dynamics model with; transition matrix A,
control matrix B, and symmetric zero-mean process noise with
symmetric covariance matrix W
"""
struct LinearDynamicsModel{a<:Number, b<:Number,
                c<:Number} <: DynamicsModel
    A::AbstractMatrix{a}
    B::AbstractMatrix{b}
    W::Symmetric{c}
    #=
    function LinearDynamicsModel{a,b,c}(A::Matrix{a},B::Matrix{b},W::Symmetrix{c})
        where{a<:Number,b<:Number,c<:Number}

        (size(A,1) == size(B,1) == size(W,1)) ?
        new(A,B,W) : error("first dimensions not matching")
    end
    =#
end

function LinearDynamicsModel(A::AbstractMatrix, B::AbstractMatrix,
                             W::AbstractMatrix)

    return LinearDynamicsModel(A,B,Symmetric(W))
end

"""
    LinearObservationModel(C::AbstractMatrix,D::AbstractMatrix,V::Symmetric)
    LinearObservationModel(C::AbstractMatrix,D::AbstractMatrix,V::AbstractMatrix)
    LinearObservationModel(C::AbstractMatrix,V::Symmetric)
    LinearObservationModel(C::AbstractMatrix,V::AbstractMatrix)

Construct linear observation dynamics model with; transition matrix C,
control matrix B, and symmetric zero-mean measurement noise with
symmetric covariance matrix V
"""
struct LinearObservationModel{a<:Number,b<:Number,c<:Number} <: ObservationModel
    C::AbstractMatrix{a}
    D::AbstractMatrix{b}
    V::Symmetric{c}
    #=
    LinearObservationModel(C,D,V) = (size(C,1) == size(D,1) == size(V,1)) ?
        new(C,D,V) : error("first dimensions not matching")
    =#
end

function LinearObservationModel(C::AbstractMatrix, D::AbstractMatrix,
                                V::AbstractMatrix)
    return LinearObservationModel(C,D,Symmetric(V))
end

function LinearObservationModel(C::AbstractMatrix, V::Symmetric)
    n = size(C,1)
    return LinearObservationModel(C,zeros(Bool,n,n),V)
end

function LinearObservationModel(C::AbstractMatrix, V::AbstractMatrix)
    n = size(C,1)
    return LinearObservationModel(C,zeros(Bool,n,n),Symmetric(V))
end

"""
    NonlinearDynamicsModel(f::Function,W::Symmetric)
    NonlinearDynamicsModel(f::Function,W::AbstractMatrix)

Construct nonlinear dynamics model with transition function f
and symmetric zero-mean process noise with symmetric covariance matrix W
"""
struct NonlinearDynamicsModel{c<:Number} <: DynamicsModel
    f::Function
    W::Symmetric{c}
end

function NonlinearDynamicsModel(f::Function, W::AbstractMatrix)
    return NonlinearDynamicsModel(f,Symmetric(W))
end

"""
    NonlinearObservationModel(h::Function,V::Symmetric)
    NonlinearObservationModel(h::Function,V::AbstractMatrix)

Construct nonlinear observation dynamics model with measurement function h
and symmetric zero-mean measurement noise with symmetric covariance matrix V
"""
struct NonlinearObservationModel{c<:Number} <: ObservationModel
    h::Function
    V::Symmetric{c}
end

function NonlinearObservationModel(h::Function, V::AbstractMatrix)
    return NonlinearObservationModel(h,Symmetric(V))
end

### Filters ###

"""
`AbstractFilter` is an abstract type to encapsulate different kinds
of discrete gaussian filters
"""
abstract type AbstractFilter end


"""
    KalmanFilter(d::LinearDynamicsModel,o::LinearObservationModel)

Construct Kalman filter with LinearDynamicsModel d and
LinearObservationModel o.
"""
struct KalmanFilter <: AbstractFilter
    d::LinearDynamicsModel
    o::LinearObservationModel
end

"""
    ExtendedKalmanFilter(d::DynamicsModel,o::ObservationModel)
    KalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number,
        α::Float,β::Float)

Construct Extended Kalman filter with DynamicsModel d and
ObservationModel o.
"""
struct ExtendedKalmanFilter <: AbstractFilter
    d::DynamicsModel
    o::ObservationModel
    ExtendedKalmanFilter(d,o) =
        (d isa LinearDynamicsModel && o isa LinearObservationModel) ?
        error("linear models: use kf over ekf") : new(d,o)
end

"""
    UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number,
        α::Float,β::Float)
    UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::Number)
    UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel)

Construct Unscented Kalman filter with DynamicsModel d, ObservationModel o,
and UKF parameters λ, α, and β. Default constructor uses α/β formulation from
Probabilistic Robotics, second constructor reduces complexity, third
constructor defaults λ to 2, as is commonly done.
"""
struct UnscentedKalmanFilter{a<:Number,b<:Number,c<:Number} <: AbstractFilter
    d::DynamicsModel
    o::ObservationModel
    λ::a
    α::b
    β::c
    #=
    UnscentedKalmanFilter(d,o,λ,α,β) = (d isa LinearDynamicsModel && o isa LinearObservationModel) ?
        error("linear models: use kf over ukf") : new(d,o,λ,α,β)
    =#
end

function UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::a) where a<:Number
    return UnscentedKalmanFilter{a,Int8,Int8}(d,o,λ,1,0)
end

function UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel)
    return UnscentedKalmanFilter{Int8,Int8,Int8}(d,o,2,1,0)
end

### Belief States ###

"""
    GaussianBelief(μ::AbstractVector,Σ::Symmetric)
    GaussianBelief(μ::AbstractVector,Σ::AbstractMatrix)

Construct a gaussian belief, consisting of mean vector μ
and symmetric covariance matrix Σ
"""
struct GaussianBelief{a<:Number,b<:Number}
    μ::AbstractVector{a}
    Σ::Symmetric{b}
    #=
    function GaussianBelief{a,b}(μ::Vector{a},
        Σ::Symmetric{b}) where {a<:Number,b<:Number}
        (size(μ,1) == size(Σ,1)) ?
        new(μ,Σ) : error("first dimensions not matching")
    end
    =#
end

function GaussianBelief(μ::AbstractVector,Σ::AbstractMatrix)
    return GaussianBelief(μ,Symmetric(Σ))
end
