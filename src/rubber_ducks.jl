__precompile__(true)
module RubberDucks

export rubber_ducks_in_earth
"""
Computes the number of rubber ducks that would fill the volume of the Earth.

Assumes a spherical Earth, and liquid rubber ducks (to avoid the packing problem).

Arguments:
- `num_ruber::Real` is the first function argument. [_m^3_]
- `r_earth::Real` is the first function argument. Default: 6378000 [_m_]

Returns:
- `num_ducks::Float64` Number of rubber ducks that could 
"""
function rubber_ducks_in_earth(duck_volume::Real, r_earth::Real=6378000)
    # Volume of spherrical Earth
    v_earth =  4.0/3.0*pi*r_earth^3

    # Return number of ducks in Earth
    return v_earth/duck_volume
end

export RUBBER_DUCK_PER_SEC
"""
Number of rubber ducks produced per second in the entire world on average.

Units: _ducks/sec_
"""
const RUBBER_DUCK_PER_SEC = 200.0

export add_one
"""
Adds one to a value

Argments:
- `val:Real` A value

Returns:
- `x::Float64` Value plus one
"""
function add_one(val::Real)
    return val + 1
end

end # End of module RubberDucks (This comment is nice to distinguish from an end)