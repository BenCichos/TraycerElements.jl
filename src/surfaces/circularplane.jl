@kwdef struct CircularPlane{N} <: Surface{N}
    normal::SVector{N,Float64}
    center::SVector{N,Float64}
    radius::Float64
    interface::Interface

    function CircularPlane(center::SVector{2,<:Real}, normal::SVector{2,<:Real}, radius::Real, interface)
        new{2}(normalize(normal), center, radius, interface)
    end

    function CircularPlane(center::SVector{3,<:Real}, normal::SVector{3,<:Real}, radius::Real, interface)
        new{3}(normalize(normal), center, radius, interface)
    end
end

normal(circularplane::CircularPlane) = circularplane.normal
center(circularplane::CircularPlane) = circularplane.center
radius(circularplane::CircularPlane) = circularplane.radius

show(io::IO, circularplane::CircularPlane) = print(io, "CircularPlane($(center(circularplane)), $(normal(circularplane)), $(radius(circularplane)))")

@approx function minintersection!(minintersection::MinIntersection, circularplane::CircularPlane{N}, ray::Ray{N}) where {N}
    plane_center = center(circularplane)
    plane_normal = normal(circularplane)

    nominator = dot(plane_center - origin(ray), plane_normal)
    iszero(nominator) && return nothing

    denominator = dot(plane_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator
    0 < alpha || return nothing
    point_of_intersection = origin(ray, alpha)
    axis_aligned = invquaternionz(plane_normal) * (point_of_intersection - plane_center)
    norm(axis_aligned) <= radius(circularplane) && return distance!(minintersection, alpha)
    return nothing
end
