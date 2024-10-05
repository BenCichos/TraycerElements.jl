struct Bubble{N} <: CompoundObject{N}
    primitives::Tuple{Sphere{N},Sphere{N}}

    function Bubble(center::SVector{N,<:Real}, outerradius::Real, innerradius::Real, outer_interface, inner_interface) where {N}
        outer_sphere = Sphere(center, outerradius, outer_interface)
        inner_sphere = Sphere(center, innerradius, inner_interface)
        new{N}(outer_sphere, inner_sphere)
    end
end
