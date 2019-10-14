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
struct LinearDynamicsModel{a<:AbstractMatrix, b<:AbstractMatrix,
                c<:Symmetric} <: DynamicsModel
    A::a
    B::b
    W::c
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
    predict(m::LinearDynamicsModel, x::AbstractVector{<:Number}, u::AbstractVector{<:Number})
    predict(m::LinearDynamicsModel, x::AbstractVector{<:Number}, u::AbstractVector{<:Number}, rng::AbstractRNG)

Uses the linear dynamics model to propagate the state x one step forward in time with control input u.
If rng is given, it adds process noise. 
"""
function predict(m::LinearDynamicsModel, x::AbstractVector{<:Number}, 
                 u::AbstractVector{<:Number})
    return m.A * x + m.B * u
end
function predict(m::LinearDynamicsModel, x::AbstractVector{T}, 
                 u::AbstractVector{T}, 
                 rng::AbstractRNG) where T<:Number
    return m.A * x + m.B * u + cholesky(m.W).L * randn(rng, T, size(m.W, 1))
end

function jacobian(m::LinearDynamicsModel, x, u)
    return m.A
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
struct LinearObservationModel{a<:AbstractMatrix,b<:AbstractMatrix,c<:Symmetric} <: ObservationModel
    C::a
    D::b
    V::c
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
    measure(m::LinearObservationModel, x::AbstractVector{<:Number}, u::AbstractVector{<:Number})
    measure(m::LinearObservationModel, x::AbstractVector{T}, u::AbstractVector{T}, rng::AbstractRNG) where T<:Number

Returns an observation of state x according to the linear observation model m, with control inputs u.
If rng is passed, adds additive Gaussian noise to the observation.
"""
function measure(m::LinearObservationModel, x::AbstractVector{<:Number}, 
                 u::AbstractVector{<:Number})
    return m.C * x + m.D * u
end
function measure(m::LinearObservationModel, x::AbstractVector{T}, 
                 u::AbstractVector{T}, 
                 rng::AbstractRNG) where T<:Number
    return m.C * x + m.D * u + cholesky(m.V).L * randn(rng, T, size(m.V, 1))
end

function jacobian(m::LinearObservationModel, x, u)
    m.C
end

"""
    NonlinearDynamicsModel(f::Function,W::Symmetric)
    NonlinearDynamicsModel(f::Function,W::AbstractMatrix)

Construct nonlinear dynamics model with transition function f
and symmetric zero-mean process noise with symmetric covariance matrix W
"""
struct NonlinearDynamicsModel{F<:Function, c<:Symmetric} <: DynamicsModel
    f::F
    W::c
end

function NonlinearDynamicsModel(f::Function, W::AbstractMatrix)
    return NonlinearDynamicsModel(f,Symmetric(W))
end

"""
    predict(m::NonLinearDynamicsModel, x::AbstractVector{<:Number}, u::AbstractVector{<:Number})
    predict(m::NonLinearDynamicsModel, x::AbstractVector{<:Number}, u::AbstractVector{<:Number}, rng::AbstractRNG)

Uses the non linear dynamics model to propagate the state x one step forward in time with control input u.
If rng is given, it adds process noise. 
"""
function predict(m::NonlinearDynamicsModel, x::AbstractVector{<:Number}, 
                 u::AbstractVector{<:Number})
    return m.f(x, u)
end
function predict(m::NonlinearDynamicsModel, x::AbstractVector{T}, 
                 u::AbstractVector{T}, 
                 rng::AbstractRNG) where T<:Number
    return m.f(x, u) + cholesky(m.W).L * randn(rng, T, size(m.W, 1))
end

function jacobian(m::NonlinearDynamicsModel, x::AbstractVector, u::AbstractVector)
    return ForwardDiff.jacobian(μ -> m.f(μ, u), x)
end

"""
    NonlinearObservationModel(h::Function,V::Symmetric)
    NonlinearObservationModel(h::Function,V::AbstractMatrix)

Construct nonlinear observation dynamics model with measurement function h
and symmetric zero-mean measurement noise with symmetric covariance matrix V
"""
struct NonlinearObservationModel{F<:Function, c<:Symmetric} <: ObservationModel
    h::F
    V::c
end

function NonlinearObservationModel(h::Function, V::AbstractMatrix)
    return NonlinearObservationModel(h,Symmetric(V))
end

"""
    measure(m::LinearObservationModel, x::AbstractVector{<:Number}, u::AbstractVector{<:Number})
    measure(m::LinearObservationModel, x::AbstractVector{T}, u::AbstractVector{T}, rng::AbstractRNG) where T<:Number

Returns an observation of state x according to the non linear observation model m, with control inputs u.
If rng is passed, adds additive Gaussian noise to the observation.
"""
function measure(m::NonlinearObservationModel, x::AbstractVector{<:Number}, 
                 u::AbstractVector{<:Number})
    return m.h(x, u)
end
function measure(m::NonlinearObservationModel, x::AbstractVector{T}, 
                 u::AbstractVector{T}, 
                 rng::AbstractRNG) where T<:Number
    return m.h(x, u) + cholesky(m.V).L * randn(rng, T, size(m.V, 1))
end

function jacobian(m::NonlinearObservationModel, x::AbstractVector, u::AbstractVector)
    return ForwardDiff.jacobian(μ -> m.h(μ, u), x)
end
