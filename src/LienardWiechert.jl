module LienardWiechert

using ArgCheck
using LinearAlgebra

export SpaceTimeCurve


"Discrete space-time curve with times t and coordinatex x,y,z."
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

end
