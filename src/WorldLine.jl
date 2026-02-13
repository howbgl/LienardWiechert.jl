
export WorldLine


"Discrete space-time curve with times t and space coordinates x,y, and z."


"""    WorldLine{T}(data::Matrix{T})

Discrete space-time curve with times t and space coordinates x,y, and z. 

# Fields
- `data::Matrix{T}`: A 7xN matrix holding with the columns [t, x, y, z, dx/dt, dy/dt, dz/dt]

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
struct WorldLine{T}
    data::Matrix{T}
    function WorldLine{T}(data::Matrix{T}) where T
        @argcheck size(data, 2) == 7 "Data matrix must have 7 columns (t, x, y, z, dx/dt, dy/dt, dz/dt)."
        @argcheck all(diff(data[:, 1]) .> 0) "Time vector (first column) must have strictly increasing values."
        @argcheck all(sum(data[:, 5:7].^2, dims=2) .< 1) "Speed must be less than the speed of light (c=1)."
        return new{T}(data)
    end
end

function WorldLine(data::Matrix{<:Real})
    if size(data, 2) == 4
        _data = Matrix{eltype(data)}(undef, size(data, 1), 7)
        _data[:, 1:4] .= data
        speed_finite_diff!(_data)
        return WorldLine{eltype(data)}(_data)
    else # Internal constructor catches wrong column size, just call it here.
        return WorldLine{eltype(data)}(data)
    end
end

function WorldLine(
    t::Vector{T},
    x::Vector{T},
    y::Vector{T},
    z::Vector{T},
    dxdt::Vector{T},
    dydt::Vector{T},
    dzdt::Vector{T}) where T

    @argcheck length(t) == length(x) == length(y) == length(z) == length(dxdt) == length(dydt) == length(dzdt) "All input vectors must have the same length."

    data = Matrix{T}(undef, length(t), 7)
    for (i, col) in enumerate((t, x, y, z, dxdt, dydt, dzdt))
        data[:, i] .= col
    end
    return WorldLine{T}(data)
end

function WorldLine(t::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}) where T
    @argcheck length(t) == length(x) == length(y) == length(z) "All input vectors must have the same length."
    data = Matrix{T}(undef, length(t), 7)
    data[:, 1] = t
    data[:, 2] = x
    data[:, 3] = y
    data[:, 4] = z
    speed_finite_diff!(data)
    return WorldLine{T}(data)
end

WorldLine(args::Vararg{AbstractVector{<:Real}}) = WorldLine(promote(collect.(args)...)...)

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

function speed_finite_diff!(data::Matrix{<:Real})
    @warn "No speed data provided. Falling back to finite difference approximation."

    dt = diff(data[:, 1]) # Time differences
    data[1:end-1, 5:7] .= diff(data[:, 2:4], dims=1) # Compute spatial diffs
    for i in 1:3
        data[1:end-1, i+4] ./= dt # Divide spatial diffs by time diffs to get speed
    end
    data[end, 5:7] .= data[end-1, 5:7] # Set last speed equal to second last (maybe NaN?)
    return nothing
end

function _below_lightspeed(points::Matrix{<:Real})
    diffs = diff(points, dims=1)
    norms = sqrt.(sum(diffs[:, 2:end].^2, dims=2))
    return all(norms .< diffs[:, 1]) # Check if speed < c (=1)
end
