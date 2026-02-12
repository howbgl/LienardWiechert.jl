
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
        @argcheck _below_lightspeed(data) "WorldLine must be below the speed of light (c=1)."
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

WorldLine(args::Vararg{AbstractVector{<:Real}, 4}) = WorldLine(promote(collect.(args)...)...)

##### Base methods

Base.:(==)(c1::WorldLine, c2::WorldLine) = c1.data == c2.data

Base.getproperty(c::WorldLine, sym::Symbol) = 
    if sym == :t
        c.data[:, 1]
    elseif sym == :x
        c.data[:, 2]
    elseif sym == :y
        c.data[:, 3]
    elseif sym == :z
        c.data[:, 4]
    else
        getfield(c, sym)
    end

function Base.isapprox(
    wl1::WorldLine{T},
    wl2::WorldLine{U};
    atol::Real=0,
    rtol=atol>0 ? 0 : √eps(promote_type(T,U)),
    nans::Bool=false) where {T,U}
    
    @argcheck size(wl1.data) == size(wl2.data) "WorldLines must have the same length."
    return all(isapprox.(wl1.data, wl2.data; atol=atol, rtol=rtol, nans=nans))
end

function Base.show(io::IO, wl::WorldLine{T}) where T
    print(io, "WorldLine{$T} with $(size(wl.data, 1)) points.")
end

function Base.show(io::IO, ::MIME"text/plain", wl::WorldLine{T}) where T
    println(io, 
    "WorldLine{$T} with $(size(wl.data, 1)) points. Columns are [t, x, y, z] with data:")

    recur_io = IOContext(io, :SHOWN_SET => wl)
    Base.print_array(recur_io, wl.data)
end

### Helper functions

function _below_lightspeed(points::Matrix{<:Real})
    diffs = diff(points, dims=1)
    norms = sqrt.(sum(diffs[:, 2:end].^2, dims=2))
    return all(norms .< diffs[:, 1]) # Check if speed < c (=1)
end
