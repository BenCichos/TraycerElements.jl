@kwdef struct LinkedElement{N,PE<:PrimitiveElement{N},LE<:PrimitiveElement{N}} <: CompoundElement{N}
    primitive_element::PE
    linked_element::LE
end
export LinkedElement

primitive_element(le::LinkedElement) = le.primitive_element
linked_element(le::LinkedElement) = le.linked_element
tracingtypes(linked_element::LinkedElement) = [linked_element,]
primitives(linked_element::LinkedElement) = [linked_element.primitive_element, linked_element.linked_element]
interface(linked_element::LinkedElement) = interface(primitive_element(linked_element))
isobject(linked_element::LinkedElement) = isobject(primitive_element(linked_element))
onintersect(linked_element::LinkedElement{N}, ray::Ray{N}, distance::Float64, normal::SVector{N,Float64}) where {N} = onintersect(linked_element.linked_element, ray, distance, normal)

normal(linked_element::LinkedElement{N}, point::SVector{N,Float64}) where {N} = normal(primitive_element(linked_element), point)
onsurface(linked_element::LinkedElement, point::SVector{N,Float64}) where {N} = onsurface(primitive_element(linked_element), point)
minintersection!(minintersection::MinIntersection, linked_element::LinkedElement{N}, ray::Ray{N}) where {N} = minintersection!(minintersection, primitive_element(linked_element), ray)
