@kwdef struct Ring{N} <: Surface{N}
    center::SVector{N,Float64}
    normal::SVector{N,Float64}
    innerradius::Float64
    outerradius::Float64
    interface::Interface
    opening_angle_cos::Float64
    opening_rotation::Float64
    opening_rotation_vec::SVector{N,Float64}

    function Ring(center::SVector{3,<:Real}, normal::SVector{3,<:Real}, innerradius::Real, outerradius::Real, interface::Interface; opening_angle::Real=pi, opening_rotation::Real=0)
        new{3}(center, normal, innerradius, outerradius, interface, cos(opening_angle), opening_rotation, quaternion(AXIS3_Z, opening_rotation) * AXIS3_X)
    end
end

center(r::Ring) = r.center
normal(r::Ring) = r.normal
innerradius(r::Ring) = r.innerradius
outerradius(r::Ring) = r.outerradius
opening_angle(r::Ring) = acos(r.opening_angle_cos)
opening_angle_cos(r::Ring) = r.opening_angle_cos
opening_rotation(r::Ring) = r.opening_rotation
opening_rotation_vec(r::Ring) = r.opening_rotation_vec

normal(r::Ring{N}, ::SVector{N,Float64}) where {N} = normal(r)

@approx function minintersection!(minintersection::MinIntersection, ring::Ring{N}, ray::Ray{N}) where {N}
    ring_center = center(ring)
    ring_normal = normal(ring)

    nominator = dot(ring_center - origin(ray), ring_normal)
    iszero(nominator) && return nothing

    denominator = dot(ring_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator
    0 < alpha || return nothing
    point_of_intersection = origin(ray, alpha)
    axis_aligned = invquaternionz(ring_normal) * (point_of_intersection - ring_center)
    angle = dot(normalize(axis_aligned), opening_rotation_vec(ring))
    (angle <= abs(opening_angle_cos(ring))) && return nothing


    (innerradius(ring) <= norm(axis_aligned) <= outerradius(ring)) && return distance!(minintersection, alpha)
    return nothing
end
