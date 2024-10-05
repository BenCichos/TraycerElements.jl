@kwdef struct Sphere{N} <: Object{N}
    center::SVector{N,Float64}
    radius::Float64
    interface::Interface
end


radius(sphere::Sphere) = sphere.radius
center(sphere::Sphere) = sphere.center
normal(sphere::Sphere{N}, point::SVector{N,Float64}) where {N} = normalize(point - center(sphere))

show(io::IO, sphere::Sphere) = print(io, "Sphere($(center(sphere)), $(radius(sphere)))")

@approx onsurface(sphere::Sphere{N}, point::SVector{N}) where {N} = norm(point - center(sphere)) == radius(sphere)

@approx function minintersection!(minintersection::MinIntersection, sphere::Sphere{N}, ray::Ray{N}) where {N}
    distance_apart = origin(ray) - center(sphere)
    ray_direction = direction(ray)
    sphere_radius = radius(sphere)

    a = dot(ray_direction, ray_direction)
    b = 2 * dot(ray_direction, distance_apart)
    c = dot(distance_apart, distance_apart) - sphere_radius^2

    alpha_1, alpha_2 = findroots(a, b, c)

    0 < alpha_2 || return nothing
    0 < alpha_1 || return distance!(minintersection, alpha_2)
    distance!(minintersection, alpha_1)
end
