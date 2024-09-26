struct Sphere{N} <: Object{N}
    center::SVector{N,Float64}
    radius::Float64
    interface::Interface

    function Sphere(center::SVector{2,<:Real}, radius::Real, interface)
        new{2}(center, radius, interface)
    end

    function Sphere(center::SVector{3,<:Real}, radius::Real, interface)
        new{3}(center, radius, interface)
    end
end

function Sphere(center::Vector{<:Real}, radius::Real, refractive_index::Real)
    Sphere(SVector{length(center)}(center), radius, Interface(refractive_index))
end

function Sphere(; center::SVector{N,<:Real}, radius::Real, refractive_index::Real) where {N}
    Sphere(center, radius, Interface(refractive_index))
end

radius(sphere::Sphere) = sphere.radius
center(sphere::Sphere) = sphere.center
normal(sphere::Sphere{N}, point::SVector{N,Float64}) where {N} = normalize(point - center(sphere))

show(io::IO, sphere::Sphere) = print(io, "Sphere($(center(sphere)), $(radius(sphere)))")

function minintersection!(minintersection::MinIntersection, sphere::Sphere{N}, ray::Ray{N}) where {N}
    distance_apart = origin(ray) - center(sphere)
    ray_direction = direction(ray)
    sphere_radius = radius(sphere)

    a = dot(ray_direction, ray_direction)
    b = 2 * dot(ray_direction, distance_apart)
    c = dot(distance_apart, distance_apart) - sphere_radius^2

    alpha_1, alpha_2 = roots(a, b, c)

    alpha_2 < 1e-10 && return nothing
    alpha_1 < 1e-10 && return distance!(minintersection, alpha_2)
    distance!(minintersection, alpha_1)
end
