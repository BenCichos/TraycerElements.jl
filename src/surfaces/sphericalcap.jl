@kwdef struct SphericalCap{N} <: Surface{N}
    center::SVector{N,Float64}
    axis::SVector{N,Float64}
    radius::Float64
    theta::Float64
    interface::Interface

    function SphericalCap(center::SVector{2,<:Real}, axis::SVector{2,<:Real}, radius::Real, theta::Real, interface)
        new{2}(center, normalize(axis), radius, theta, interface)
    end

    function SphericalCap(center::SVector{3,<:Real}, axis::SVector{3,<:Real}, radius::Real, theta::Real, interface)
        new{3}(center, normalize(axis), radius, theta, interface)
    end

end

radius(spherical_cap::SphericalCap) = spherical_cap.radius
theta(spherical_cap::SphericalCap) = spherical_cap.theta
center(spherical_cap::SphericalCap) = spherical_cap.center
axis(spherical_cap::SphericalCap) = spherical_cap.axis
normal(spherical_cap::SphericalCap{N}, point::SVector{N,Float64}) where {N} = normalize(point - center(spherical_cap))

function onsurface(spherical_cap::SphericalCap{N}, point::SVector{N,Float64}) where {N}
    cosine_normal_axis = dot(normal(spherical_cap, point), axis(spherical_cap))
    theta(spherical_cap) < acos(cosine_normal_axis) && return false
    return true
end

show(io::IO, spherical_cap::SphericalCap) = print(io, "SphericalCap($(center(spherical_cap)), $(axis(spherical_cap)), $(radius(spherical_cap)), $(theta(spherical_cap)))")

@approx function minintersection!(minintersection::MinIntersection, spherical_cap::SphericalCap{N}, ray::Ray{N}) where {N}
    distance_apart = origin(ray) - center(spherical_cap)
    ray_direction = direction(ray)
    spherical_cap_radius = radius(spherical_cap)

    a = dot(ray_direction, ray_direction)
    b = 2 * dot(ray_direction, distance_apart)
    c = dot(distance_apart, distance_apart) - spherical_cap_radius^2

    alpha_1, alpha_2 = findroots(a, b, c)

    0 < alpha_2 || return nothing

    alpha_2 = onsurface(spherical_cap, origin(ray, alpha_2)) ? alpha_2 : Inf

    0 < alpha_1 || return distance!(minintersection, alpha_2)

    onsurface(spherical_cap, origin(ray, alpha_1)) && return distance!(minintersection, alpha_1)

    distance!(minintersection, alpha_2)
end
