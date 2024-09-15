struct LightCone{N} <: Surface{N}
    ring::Ring{N}
    tan_angle::Float64
    focal_length::Float64
    intensity::Float64


    function LightCone(ring::Ring{N}, focal_length::Real, intensity::Real) where {N}
        new{N}(ring, outerradius(ring) / focal_length, focal_length, intensity)
    end

    function LightCone(center::SVector{N,<:Real}, normal::SVector{N,<:Real}, innerradius::Real, outerradius::Real, focal_length::Real; interface::Interface=ConsumingInterface(), intensity::Real=1.0) where {N}
        ring = Ring(center, normal, innerradius, outerradius, interface)
        new{N}(ring, outerradius / focal_length, focal_length, intensity)
    end
end

ring(cone::LightCone) = cone.ring
focal_length(cone::LightCone) = cone.focal_length
tan_angle(cone::LightCone) = cone.tan_angle
intensity(cone::LightCone) = cone.intensity
function normal(cone::LightCone, point::SVector{N,Float64}) where {N}
    delta_focal_length = (outerradius(ring(cone)) - norm(point - center(ring(cone)))) / tan_angle(cone)
    focal_point = center(ring(cone)) + (focal_length(cone) - delta_focal_length) * normal(ring(cone))
    return normalize(focal_point - point)
end
interface(cone::LightCone) = interface(ring(cone))
minintersection!(minintersection::MinIntersection, lightplane::LightCone{N}, ray::Ray{N}) where {N} = minintersection!(minintersection, ring(lightplane), ray)

function pointnormals(cone::LightCone; angularresolution::Int=48, radialresolution::Int=10)
    points = [SA[0, r*sin(theta), r*cos(theta)] + center(ring(cone)) for theta in LinRange(0, 2pi, angularresolution) for r in LinRange(innerradius(ring(cone)), outerradius(ring(cone)), radialresolution)]
    normals = [normal(cone, point) for point in points]
    return (points, normals)
end
export pointnormals
