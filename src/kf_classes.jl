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
    LinearDynamicsModel(A::Matrix,B::Matrix,W::Symmetric)
    LinearDynamicsModel(A::Matrix,B::Matrix,W::Matrix)

Construct linear dynamics model with; transition matrix A,
control matrix B, and symmetric zero-mean process noise with
symmetric covariance matrix W
"""
mutable struct LinearDynamicsModel{a,b,c} <: DynamicsModel
    A::Matrix{a}
    B::Matrix{b}
    W::Symmetric{c}
    LinearDynamicsModel(A,B,W) = (size(A,1) == size(B,1) == size(W,1)) ?
        new(A,B,W) : error("first dimensions not matching")
end

function LinearDynamicsModel(A::Matrix, B::Matrix, W::Matrix)
    return LinearDynamicsModel(A,B,Symmetric(W))
end

"""
    LinearObservationModel(C::Matrix,D::Matrix,V::Symmetric)
    LinearObservationModel(C::Matrix,D::Matrix,V::Matrix)
    LinearObservationModel(C::Matrix,V::Symmetric)
    LinearObservationModel(C::Matrix,V::Matrix)

Construct linear observation dynamics model with; transition matrix C,
control matrix B, and symmetric zero-mean measurement noise with
symmetric covariance matrix V
"""
mutable struct LinearObservationModel{a,b,c} <: ObservationModel
    C::Matrix{a}
    D::Matrix{b}
    V::Symmetric{c}
    LinearObservationModel(C,D,V) = (size(C,1) == size(D,1) == size(V,1)) ?
        new(C,D,V) : error("first dimensions not matching")
end

function LinearObservationModel(C::Matrix, D::Matrix, V::Matrix)
    return LinearObservationModel(C,D,Symmetric(V))
end

function LinearObservationModel(C::Matrix, V::Symmetric)
    n = size(C,1)
    return LinearObservationModel(C,zeros(Bool,n,n),V)
end

function LinearObservationModel(C::Matrix, D::Matrix, V::Matrix)
    n = size(C,1)
    return LinearObservationModel(C,zeros(Bool,n,n),Symmetric(V))
end

"""
    NonlinearDynamicsModel(f::Function,W::Symmetric)
    NonlinearDynamicsModel(f::Function,W::Matrix)

Construct nonlinear dynamics model with transition function f
and symmetric zero-mean process noise with symmetric covariance matrix W
"""
mutable struct NonlinearDynamicsModel{c} <: DynamicsModel
    f::Function
    W::Symmetric{c}
end

function NonlinearDynamicsModel(f::Function, W::Matrix)
    return NonlinearDynamicsModel(f,Symmetric(W))
end

"""
    NonlinearObservationModel(h::Function,V::Symmetric)
    NonlinearObservationModel(h::Function,V::Matrix)

Construct nonlinear observation dynamics model with measurement function h
and symmetric zero-mean measurement noise with symmetric covariance matrix V
"""
mutable struct NonlinearObservationModel{c} <: ObservationModel
    h::Function
    V::Symmetric{c}
end

function NonlinearObservationModel(h::Function, V::Matrix)
    return NonlinearObservationModel(h,Symmetric(V))
end

"""
    KalmanFilter(d::LinearDynamicsModel,o::LinearObservationModel)

Construct Kalman filter with LinearDynamicsModel d and
LinearObservationModel o.
"""
mutable struct KalmanFilter
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
mutable struct ExtendedKalmanFilter
    d::DynamicsModel
    o::ObservationModel
    ExtendedKalmanFilter(d,o) = (d isa LinearDynamicsModel and o isa LinearObservationModel) ?
        error("linear models: use kf over ekf"): new(d,o)
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
mutable struct UnscentedKalmanFilter{a,b,c}
    d::DynamicsModel
    o::ObservationModel
    λ::a
    α::b
    β::c
    UnscentedKalmanFilter(d,o,λ,α,β) = (d isa LinearDynamicsModel and o isa LinearObservationModel) ?
        error("linear models: use kf over ukf"): new(d,o,λ,α,β)
end

function UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel,λ::a) where a
    return UnscentedKalmanFilter{a,Int8,Int8}(d,o,λ,1,0)
end

function UnscentedKalmanFilter(d::DynamicsModel,o::ObservationModel)
    return UnscentedKalmanFilter{Int8,Int8,Int8}(d,o,2,1,0)
end
