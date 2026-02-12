
export WorldLine


"Discrete space-time curve with times t and space coordinates x,y, and z."


"""    WorldLine{T}(data::Matrix{T})

Discrete space-time curve with times t and space coordinates x,y, and z. 

# Fields
- `data::Matrix{T}`: A 4xN matrix holding the data [t, x, y, z] in rows.
"""
struct WorldLine{T}
    data::Matrix{T}
    function WorldLine{T}(data::Matrix{T}) where T
        @argcheck size(data, 2) == 4 "Data matrix must have 4 columns (t, x, y, z)."
        @argcheck all(diff(data[:, 1]) .> 0) "Time vector (first column) must have strictly increasing values."
        return new{T}(data)
    end
end

WorldLine(data::Matrix{<:Real}) = WorldLine{eltype(data)}(data)

function WorldLine(t::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}) where T
    @argcheck length(t) == length(x) == length(y) == length(z) "All input vectors must have the same length."
    data = Matrix{T}(undef, length(t), 4)
    data[:, 1] = t
    data[:, 2] = x
    data[:, 3] = y
    data[:, 4] = z
    return WorldLine{T}(data)
end

function WorldLine(args::Vararg{AbstractVector{<:Real}, 4})
    args = promote(args...)
    WorldLine{eltype(args[1])}(args...)
end

Base.:(==)(c1::WorldLine, c2::WorldLine) = c1.data == c2.data
