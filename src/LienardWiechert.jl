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
        @argcheck all(diff(t) .> 0) "Time vector must have strictly increasing values."
        new{T}(t, x, y, z)
    end
end
function SpaceTimeCurve(t::Vector{<:Real}, x::Vector{<:Real}, y::Vector{<:Real}, z::Vector{<:Real}) 
    args = promote(t, x, y, z)
    SpaceTimeCurve{eltype(args[1])}(args...)
end
function SpaceTimeCurve(args::Vararg{AbstractVector{<:Real}, 4})
    SpaceTimeCurve([collect(a) for a in args]...)
end

"Calculate the Minkowski distance squared between two spacetime points, (+---) signature."
function minkowski_distance_squared(x1,x2)
    dx = abs.(x1 .- x2)
    return minkowski_dot(dx, dx)
end

function minkowski_dot(x1, x2)
    return x1[1]*x2[1] - sum(x1[2:end] .* x2[2:end])
end

function root_linear_interpolation(x1, t1, x2, t2)
    return t1 - x1 * (t2 - t1) / (x2 - x1)
end

function retarded_time(spacetime_point, c::SpaceTimeCurve)
    @argcheck length(spacetime_point) == 4 "Spacetime point must be a 4-element vector (t, x, y, z)."
    @argcheck all(isfinite.(spacetime_point)) "Spacetime point must have finite values."
    @argcheck spacetime_point[1] > c.t[1] "Spacetime point must be in the future of the curve start."
    
    start_distance_squared = minkowski_distance_squared(
        spacetime_point, 
        [c.t[1], c.x[1], c.y[1], c.z[1]])
    
    if start_distance_squared < 0 
        throw(ArgumentError(
            "Spacetime point is outside the light cone of curve start. No retarded time solution exists."))
    end

    for i in eachindex(c.t)
        d2 = minkowski_distance_squared(spacetime_point, [c.t[i], c.x[i], c.y[i], c.z[i]])
        if d2 < 0
            # return (c.t[i-1] + c.t[i]) / 2
            d1 = minkowski_distance_squared(
                spacetime_point, 
                [c.t[i-1], c.x[i-1], c.y[i-1], c.z[i-1]])
            return root_linear_interpolation(d1, c.t[i-1], d2, c.t[i])
        end
    end


    # If we reach here, the spacetime point is outside the light cone of the entire curve
    throw(ArgumentError(
        "Spacetime point is outside the light cone of the entire curve. No retarded time solution exists."))
end

end
