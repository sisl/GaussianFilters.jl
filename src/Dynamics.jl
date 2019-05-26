### Dynamics Model ###
"""
    Dynamics(A,Q,d)
    Dynamics(A,Q)

Construct linear dynamics model with; transition matrix A,
process noise matrix Q and constant matrix d
"""
mutable struct Dynamics{A,Q,d}
A::Matrix{A}
Q::Matrix{Q}
d::Matrix{d}
end

## Constructors ##
function Dynamics(A,Q)
    return Dynamics(A,Q,d = nothing)
end

function Dynamics(A,Q,d)
    return Dynamics(A,Q,d)
end
