@kwdef struct LightSurface{N,S<:Surface{N}} <: Surface{N}
    primitive_element::S
    intensity::Float64
    focal_length::Float64
end

primitive_element(light_surface::LightSurface) = light_surface.primitive_element
interface(light_surface::LightSurface) = interface(primitive_element(light_surface))
normal(light_surface::LightSurface) = normal(primitive_element(light_surface))

minintersection!(minintersection::MinIntersection, lightplane::LightSurface{N}, ray::Ray{N}) where {N} = minintersection!(minintersection, primitive_element(light_surface), ray)
