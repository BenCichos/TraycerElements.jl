@kwdef struct Plane{N,M} <: Surface{N}
    center::SVector{N,Float64}
    normal::SVector{N,Float64}
    size::NTuple{M,Float64}
    interface::Interface

    function Plane(center::SVector{2,<:Real}, normal::SVector{2,<:Real}, size::Real=Inf, interface=NULL_INTERFACE)
        new{2,1}(center, normalize(normal), (size,), interface)
    end

    function Plane(center::SVector{3,<:Real}, normal::SVector{3,<:Real}, size::Tuple{<:Real,<:Real}=(Inf, Inf), interface=NULL_INTERFACE)
        new{3,2}(center, normalize(normal), size, interface)
    end

end

center(plane::Plane) = plane.center
normal(plane::Plane) = plane.normal
size(plane::Plane{2}) = only(plane.size)
size(plane::Plane{3}) = plane.size
hassize(plane::Plane) = !all(isinf, size(plane))

isinside(plane::Plane{2}, distance::Float64) = all(abs.(distance) .< (size(plane) ./ 2))

function minintersection!(minintersection::MinIntersection, plane::Plane{2}, ray::Ray{2})
    plane_point = center(plane)
    plane_normal = normal(plane)

    nominator = dot(plane_point - origin(ray), plane_normal)
    iszero(nominator) && return nothing

    denominator = dot(plane_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator

    alpha < 1e-10 && return nothing
    hassize(plane) || return nothing

    point_of_intersection = origin(ray, alpha)
    plane_vector = perpto(plane_normal)
    distance_from_center = dot(point_of_intersection - plane_point, plane_vector) |> abs
    distance_from_center > (size(plane) / 2) && return nothing

    distance!(minintersection, alpha)
end

function minintersection!(minintersection::MinIntersection, plane::Plane{3}, ray::Ray{3})
    plane_center = center(plane)
    plane_normal = normal(plane)

    nominator = dot(plane_center - origin(ray), plane_normal)
    iszero(nominator) && return nothing

    denominator = dot(plane_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator

    alpha < _TOLERANCE && return Inf
    hassize(plane) || return nothing

    point_of_intersection = origin(ray, alpha)
    axis_aligned_vector = quaternionz(plane_normal) * (point_of_intersection - plane_center)
    abs(axis_aligned_vector[2]) > (y_limit / 2) && return nothing

    distance!(minintersection, alpha)
end

show(io::IO, plane::Plane) = print(io, "Plane($(center(plane)), $(normal(plane)), $(size(plane))")
