struct LightPlane{N} <: Surface{N}
    plane::Plane{N}
    intensity::Float64

    function LightPlane(origin::SVector{2,<:Real}, normal::SVector{2,<:Real}, size::Real=Inf, intensity::Real=1.0, interface::Interface=ConsumingInterface())
        plane = Plane(origin, normal, interface, size)
        new{2}(plane, intensity)
    end

    function LightPlane(origin::SVector{3,<:Real}, normal::SVector{3,<:Real}, size::Tuple{<:Real,<:Real}=(Inf, Inf), intensity::Real=1.0, interface::Interface=ConsumingInterface())
        plane = Plane(origin, normal, interface, size)
        new{3}(plane, intensity)
    end
end
export LightPlane

plane(light::LightPlane) = light.plane
intensity(light::LightPlane) = light.intensity
interface(light::LightPlane) = interface(plane(light))
normal(light::LightPlane{N}) where {N} = normal(plane(light))

minintersection!(minintersection::MinIntersection, lightplane::LightPlane{N}, ray::Ray{N}) where {N} = minintersection!(minintersection, plane(lightplane), ray)
