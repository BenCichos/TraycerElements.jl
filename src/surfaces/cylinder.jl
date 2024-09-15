struct Cylinder{N} <: Surface{N}
    center::SVector{N,Float64}
    axis::SVector{N,Float64}
    height::Float64
    radius::Float64
    interface::Interface

    function Cylinder(center::SVector{2,<:Real}, axis::SVector{2,<:Real}, height::Real, radius::Real, interface)
        new{2}(center, normalize(axis), height, radius, interface)
    end

    function Cylinder(center::SVector{3,<:Real}, axis::SVector{3,<:Real}, height::Real, radius::Real, interface)
        new{3}(center, normalize(axis), height, radius, interface)
    end
end

function Cylinder(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, height::Real, radius::Real, refractive_index::Real) where {N}
    Cylinder(center, axis, height, radius, Interface(refractive_index))
end

function Cylinder(center::Vector{<:Real}, axis::Vector{<:Real}, height::Real, radius::Real, refractive_index::Real)
    Cylinder(SVector{length(center)}(center), SVector{length(axis)}(axis), height, radius, Interface(refractive_index))
end

function Cylinder(center::Vector{<:Real}, axis::Vector{<:Real}, height::Real, radius::Real, interface)
    Cylinder(SVector{length(center)}(center), SVector{length(axis)}(axis), height, radius, interface)
end

function Cylinder(; center::SVector{N,<:Real}, axis::SVector{N,<:Real}, height::Real, radius::Real, refractive_index::Float64) where {N}
    Cylinder(center, axis, height, radius, refractive_index)
end

center(cylinder::Cylinder) = cylinder.center
radius(cylinder::Cylinder) = cylinder.radius
height(cylinder::Cylinder) = cylinder.height
axis(cylinder::Cylinder) = cylinder.axis


normal(cylinder::Cylinder{N}, point::SVector{N,Float64}) where {N} = normalize(point - center(cylinder) + dot(axis(cylinder), point - center(cylinder)) * axis(cylinder))

show(io::IO, cylinder::Cylinder) = print(io, "Cylinder($(center(cylinder)), $(axis(cylinder)), $(height(cylinder)), $(radius(cylinder)))")

function minintersection!(minintersection::MinIntersection, cylinder::Cylinder{2}, ray::Ray{2})
    n = direction(ray)
    b = center(cylinder) - origin(ray)
    a = axis(cylinder)
    r = radius(cylinder)

    n_cross_a = cross(n, a)
    b_cross_a = cross(b, a)

    distance_1 = (b_cross_a - r) / n_cross_a
    distance_2 = (b_cross_a + r) / n_cross_a

    distance_2 < 1e-10 && return nothing
    alpha = distance_1 < 1e-10 ? distance_2 : distance_1
    distance_from_center = dot(a, origin(ray, alpha) - center(cylinder))
    0 <= distance_from_center <= height(cylinder) && return distance!(minintersection, alpha)

    return nothing
end

function minintersection!(minintersection::MinIntersection, cylinder::Cylinder{3}, ray::Ray{3})
    n = direction(ray)
    b = center(cylinder) - origin(ray)
    a = axis(cylinder)
    r = radius(cylinder)

    n_cross_a = cross(n, a)
    b_cross_a = cross(b, a)

    dot_n_cross_a = dot(n_cross_a, n_cross_a)

    discriminant = dot_n_cross_a * r^2 - (dot(b, n_cross_a))^2
    discriminant < 0.0 && return nothing

    sqrt_discriminant = sqrt(discriminant)

    preamble = dot(n_cross_a, b_cross_a)

    distance_1 = (preamble - sqrt_discriminant) / dot_n_cross_a
    distance_2 = (preamble + sqrt_discriminant) / dot_n_cross_a

    distance_2 < 1e-10 && return nothing
    alpha = distance_1 < 1e-10 ? distance_2 : distance_1
    distance_from_center = dot(a, origin(ray, alpha) - center(cylinder))
    0 <= distance_from_center <= height(cylinder) && return distance!(minintersection, alpha)

    return nothing
end
