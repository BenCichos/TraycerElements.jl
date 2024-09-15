struct Plane{N} <: Surface{N}
    center::SVector{N,Float64}
    normal::SVector{N,Float64}
    interface::Interface
    size::Tuple{Float64,Float64}

    function Plane(center::SVector{2,<:Real}, normal::SVector{2,<:Real}, interface, size::Real=Inf)
        new{2}(center, normalize(normal), interface, (size, 0))
    end

    function Plane(center::SVector{3,<:Real}, normal::SVector{3,<:Real}, interface, size::Tuple{<:Real,<:Real}=(Inf, Inf))
        new{3}(center, normalize(normal), interface, size)
    end

end

function Plane(center::Vector{<:Real}, normal::Vector{<:Real}, refractive_index::Real, size::Real)
    Plane(SVector{length(center)}(center), SVector{length(normal)}(normal), Interface(refractive_index), size)
end

function Plane(center::Vector{<:Real}, normal::Vector{<:Real}, interface, size::Tuple{Vararg{<:Real}})
    Plane(SVector{length(center)}(point), SVector{length(normal)}(normal), interface, size)
end

function Plane(center::SVector{2,<:Real}, normal::SVector{2,<:Real}, refractive_index::Real, size::Real=Inf)
    Plane(center, normal, Interface(refractive_index), size)
end

function Plane(center::SVector{3,<:Real}, normal::SVector{3,<:Real}, refractive_index::Real, size::Tuple{<:Real,<:Real}=(Inf, Inf))
    Plane(center, normal, Interface(refractive_index), size)
end

center(plane::Plane) = plane.center
normal(plane::Plane) = plane.normal
size(plane::Plane) = plane.size

show(io::IO, plane::Plane) = print(io, "Plane($(center(plane)), $(normal(plane)), $(size(plane))")

function minintersection!(minintersection::MinIntersection, plane::Plane{2}, ray::Ray{2})
    plane_point = center(plane)
    plane_normal = normal(plane)

    nominator = dot(plane_point - origin(ray), plane_normal)
    iszero(nominator) && return nothing

    denominator = dot(plane_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator

    alpha < 1e-10 && return nothing
    plane_size, _ = size(plane)
    isinf(plane_size) && return alpha

    distance_limit = plane_size / 2
    point_of_intersection = origin(ray, alpha)
    plane_vector = SVector(plane_normal[2], -plane_normal[1])
    distance_from_center = dot(point_of_intersection - plane_point, plane_vector) |> abs
    distance_from_center > distance_limit && return nothing

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

    alpha < 1e-10 && return Inf
    x_limit, y_limit = size(plane)
    isinf(x_limit) && isinf(y_limit) && return nothing

    rotation_matrix = rotation(plane_normal)
    point_of_intersection = origin(ray, alpha)
    axis_aligned_vector = inv(rotation_matrix) * (point_of_intersection - plane_center)
    abs(axis_aligned_vector[1]) > (x_limit / 2) && return nothing
    abs(axis_aligned_vector[2]) > (y_limit / 2) && return nothing

    distance!(minintersection, alpha)
end
