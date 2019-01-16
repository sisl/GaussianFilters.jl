__precompile__(true)
module YourSubmodule

export your_function
"""
Write something about the purpose of the function up top.

Arguments:
- `arg_1::Integer` Desribute first function argument. [input_units]
- `arg_2::Integer` Desribute second function argument. [input_units]

Returns:
- `ret_1::Float64` Describe first return value. [input_units]
- `ret_2::Float64` Describe second return value. [input_units]

Notes:
1. Make specific comments about corner cases or things to know about usage here.

References:
1. List citations to references used in writing the function here
"""
function your_function(arg_1::Integer, arg_2::Real)
    if arg_2 < 0
        return 0, 0
    end

    return arg_1 + arg_2, arg_1 - arg_2
end

end # End of module YourSubmodule (This comment is nice to distinguish from an end)