module LienardWiechert

using ArgCheck
using LinearAlgebra

export SpaceTimeCurve


"Discrete space-time curve with times t and space coordinates x,y, and z."
struct SpaceTimeCurve{T}
    t::Vector{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    function SpaceTimeCurve{T}(t::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}) where T
        @argcheck length(t) == length(x) == length(y) == length(z) "All input vectors must have the same length."
        new{T}(t, x, y, z)
    end
end
function SpaceTimeCurve(t::Vector{<:Real}, x::Vector{<:Real}, y::Vector{<:Real}, z::Vector{<:Real}) 
    args = promote(t, x, y, z)
    SpaceTimeCurve{eltype(args[1])}(args...)
end

"Calculate the Minkowski distance squared between two spacetime points, (+---) signature."
function minkowski_distance_squared(x1,x2)
    dx = x1 .- x2
    return dx[1]^2 - sum(dx[2:end].^2)
end

function find_equal_time(time::Real, c::SpaceTimeCurve)
    for i in eachindex(c.t)
        # check if time is in the interval [c.t[i], c.t[i+1]]
        if i < length(c.t) && c.t[i] <= time <= c.t[i+1]
            return i
        end
    end
    throw(ArgumentError("Time $time is out of bounds of the curve's time range."))
end

function retarded_time(spacetime_point, c::SpaceTimeCurve)
    for i in eachindex(c.t)
        r = [c.t[i], c.x[i], c.y[i], c.z[i]]
        # implement rest ...
    end
end

end
