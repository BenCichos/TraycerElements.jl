struct CircularPlane{N} <: Surface{N}
    normal::SVector{N,Float64}
    center::SVector{N,Float64}
    radius::Float64
    interface::Interface

    function CircularPlane(center::SVector{2,<:Real}, normal::SVector{2,<:Real}, radius::Real, interface)
        unit_normal = normalize(normal)
        new{2}(unit_normal, center, radius, interface)
    end

    function CircularPlane(center::SVector{3,<:Real}, normal::SVector{3,<:Real}, radius::Real, interface)
        unit_normal = normalize(normal)
        new{3}(normal, radius, dot(unit_normal, center), interface)
    end
end

function CircularPlane(center::SVector{N,<:Real}, normal::SVector{N,<:Real}, radius::Real, refractive_index::Real) where {N}
    unit_normal = normalize(normal)
    CircularPlane(center, normal, radius, Interface(refractive_index))
end

function CircularPlane(center::Vector{<:Real}, normal::Vector{<:Real}, radius::Real, refractive_index::Real)
    CircularPlane(SVector{length(center)}(center), SVector{length(normal)}(normal), radius, refractive_index)
end

function CircularPlane(center::Vector{<:Real}, normal::Vector{<:Real}, radius::Real, interface)
    CircularPlane(SVector{length(center)}(center), SVector{length(normal)}(normal), radius, interface)
end

function CircularPlane(; center::SVector{N,<:Real}, normal::SVector{N,<:Real}, radius::Real, interface) where {N}
    CircularPlane(center, normal, radius, interface)
end

normal(circularplane::CircularPlane) = circularplane.normal
center(circularplane::CircularPlane) = circularplane.center
radius(circularplane::CircularPlane) = circularplane.radius

show(io::IO, circularplane::CircularPlane) = print(io, "CircularPlane($(center(circularplane)), $(normal(circularplane)), $(radius(circularplane)))")

function minintersection!(minintersection::MinIntersection, circularplane::CircularPlane{N}, ray::Ray{N}) where {N}
    plane_center = center(circularplane)
    plane_normal = normal(circularplane)

    nominator = dot(plane_center - origin(ray), plane_normal)
    iszero(nominator) && return nothing

    denominator = dot(plane_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator
    alpha <= 1e-10 && return nothing
    point_of_intersection = origin(ray, alpha)
    axis_aligned = inv(rotation(plane_normal)) * (point_of_intersection - plane_center)
    norm(axis_aligned) <= radius(circularplane) && return distance!(minintersection, alpha)
    return nothing
end
