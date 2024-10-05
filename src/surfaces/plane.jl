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

@approx isinside(plane::Plane{2}, v::SVector{2,Float64}) = abs(v.x) < size(plane) ./ 2
@approx isinside(plane::Plane{3}, v::SVector{3,Float64}) = all(abs.(v[1:2]) .< size(plane) ./ 2)

@approx function minintersection!(minintersection::MinIntersection, plane::Plane{N}, ray::Ray{N})
    plane_center = center(plane)
    plane_normal = normal(plane)

    nominator = dot(plane_center - origin(ray), plane_normal)
    iszero(nominator) && return nothing

    denominator = dot(plane_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator

    0 < alpha || return nothing
    hassize(plane) || return nothing

    point_of_intersection = origin(ray, alpha)
    axis_aligned_vector = invquaternionz(plane_normal) * (point_of_intersection - plane_center)
    isinside(plane, axis_aligned_vector) && return nothing

    distance!(minintersection, alpha)
end

show(io::IO, plane::Plane) = print(io, "Plane($(center(plane)), $(normal(plane)), $(size(plane))")
