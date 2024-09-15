struct Ring{N} <: Surface{N}
    center::SVector{N,Float64}
    normal::SVector{N,Float64}
    innerradius::Float64
    outerradius::Float64
    interface::Interface
    opening_angle_cos::Float64
    opening_rotation::Float64
    opening_rotation_vec::SVector{N,Float64}

    function Ring(center::SVector{N,<:Real}, normal::SVector{N,<:Real}, innerradius::Real, outerradius::Real, interface::Interface; opening_angle::Real=pi, opening_rotation::Real=0) where {N}
        new{N}(center, normal, innerradius, outerradius, interface, cos(opening_angle), opening_rotation, rotate_vector(quatwithaxis(SA[0, 0, 1], opening_rotation), SA[1.0, 0.0, 0.0]))
    end
end

function quatwithaxis(axis::SVector{3,<:Real}, theta::Real)
    s, c = sincos(theta / 2)
    return Quaternion(c, s * axis[1], s * axis[2], s * axis[3])
end

function rotate_vector(q::Quaternion, u::SVector{3,<:Real})
    q_u = Quaternion(0, u[1], u[2], u[3])
    q_v = q * q_u * conj(q)
    return SVector(imag_part(q_v))
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

function minintersection!(minintersection::MinIntersection, ring::Ring{N}, ray::Ray{N}) where {N}
    ring_center = center(ring)
    ring_normal = normal(ring)

    nominator = dot(ring_center - origin(ray), ring_normal)
    iszero(nominator) && return nothing

    denominator = dot(ring_normal, direction(ray))
    iszero(denominator) && return nothing

    alpha = nominator / denominator
    alpha <= 1e-10 && return nothing
    point_of_intersection = origin(ray, alpha)
    axis_aligned = inv(rotation(ring_normal)) * (point_of_intersection - ring_center)
    angle = dot(normalize(axis_aligned), opening_rotation_vec(ring))
    (angle < abs(opening_angle_cos(ring))) && return nothing


    (innerradius(ring) <= norm(axis_aligned) <= outerradius(ring)) && return distance!(minintersection, alpha)
    return nothing
end
