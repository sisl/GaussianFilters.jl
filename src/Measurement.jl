### Measurement Model ###
"""
    Measurement(C,R)

Construct measurement model with observation matrix C and sensor
noise matrix R
"""
mutable struct Measurement(C,R)
C::Matrix{C}
R::Matrix{R}
end

## Constructor ##
function (C,R)
    return Measurement(C,R)
end
