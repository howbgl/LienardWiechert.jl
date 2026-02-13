module LienardWiechert

using ArgCheck
using LinearAlgebra

export potential, retarded_time, snapshot_2d

include("WorldLine.jl")

"Calculate the Minkowski distance squared between two spacetime points, (+---) signature."
function minkowski_distance_squared(x1,x2)
    dx = abs.(x1 .- x2)
    return minkowski_dot(dx, dx)
end

"η(r1,r2) = t1*t2 - x1*x2 - y1*y2 - z1*z2 with r1 = (t1,x1,y1,z1) and r2 = (t2,x2,y2,z2)."
function minkowski_dot(x1, x2)
    return x1[1]*x2[1] - sum(x1[2:end] .* x2[2:end])
end

"Fit a linear function to two points (x1,t1) and (x2,t2) and return the root (t when x=0)."
function root_linear_interpolation(x1, t1, x2, t2)
    return t1 - x1 * (t2 - t1) / (x2 - x1)
end

function retarded_time(spacetime_point, c::WorldLine)
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

    tr = retarded_time_nochecks(spacetime_point, c)

    if isnan(tr)
        throw(ArgumentError(
            "Curve never crosses past light cone of spacetime point. No retarded time solution exists."))        
    end
    return tr
end

function retarded_time_nochecks(spacetime_point, c::WorldLine)
    for i in eachindex(c.t)
        d2 = minkowski_distance_squared(spacetime_point, [c.t[i], c.x[i], c.y[i], c.z[i]])
        if d2 < 0
            d1 = minkowski_distance_squared(
                spacetime_point, 
                [c.t[i-1], c.x[i-1], c.y[i-1], c.z[i-1]])
            return root_linear_interpolation(d1, c.t[i-1], d2, c.t[i])
        end
    end
    return convert(eltype(c.t), NaN)
end

function interpolate_worldline(t::Real, c::WorldLine)
    @argcheck t >= c.t[1] && t <= c.t[end] "Interpolation time is outside the worldline range."
    for i in eachindex(c.t)
        if c.t[i] > t
            return c.data[i-1, :] .+ (t - c.t[i-1]) * (c.data[i, :] .- c.data[i-1, :]) / (c.t[i] - c.t[i-1])
        end
    end
end

function potential(spacetime_point, c::WorldLine)
    tr = retarded_time(spacetime_point, c)
    r  = interpolate_worldline(tr, c)
    R  = spacetime_point[2:4] .- r[2:4]
    eR = normalize(R)
    β  = r[5:7]

    denom = norm(R) * (1 - dot(eR, β))
    return [r[1] / denom; β ./ denom]
end

function snapshot_2d(xs,ys,z0,t,wl::WorldLine)
    pts = [zeros(4) for _ in xs, _ in ys]
    for i in eachindex(xs)
        println("Calculating potential for x=$(xs[i])...")
        Threads.@threads for j in eachindex(ys)
            pts[i,j] = potential([t, xs[i], ys[j], z0], wl)
        end
    end
    return pts
end

end
