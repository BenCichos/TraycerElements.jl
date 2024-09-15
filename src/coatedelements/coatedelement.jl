struct CoatedElement{N,PE<:PrimitiveElement{N},CS<:Surface{N}} <: CompoundElement{N}
    primitive_element::PE
    coating_surface::CS
    thickness::Float64
end
export CoatedElement

primitive_element(ce::CoatedElement) = ce.primitive_element
coating_surface(ce::CoatedElement) = ce.coating_surface
thickness(ce::CoatedElement) = ce.thickness

tracingtypes(coated_element::CoatedElement) = [coated_element,]
function interface(coated_element::CoatedElement, point::SVector{N,Float64}) where {N}
    onsurface(coating_surface(coated_element), point) && return CoatedInterface(interface(primitive_element(coated_element)), interface(coating_surface(coated_element)), thickness(coated_element))
    return interface(primitive_element(coated_element))
end

isobject(coated_element::CoatedElement) = isobject(primitive_element(coated_element))
primitives(coated_element::CoatedElement) = [primitive_element(coated_element), coating_surface(coated_element)]
onintersect(coated_element::CoatedElement{N}, ray::Ray{N}, distance::Float64, normal::SVector{N,Float64}) where {N} = onintersect(primitive_element(coated_element), ray, distance, normal)

normal(coated_element::CoatedElement{N}, point::SVector{N,Float64}) where {N} = normal(primitive_element(coated_element), point)
minintersection!(minintersection::MinIntersection, coated_element::CoatedElement{N}, ray::Ray{N}) where {N} = minintersection!(minintersection, primitive_element(coated_element), ray)
