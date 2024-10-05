struct Pill{N} <: CompoundObject{N}
    primitives::Tuple{SphericalCap{N},Cylinder{N},SphericalCap{N}}

    function Pill(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, radius::Real, height::Real, interface) where {N}
        front_cap = SphericalCap(center + height * axis, axis, radius, pi / 2, interface)
        cylinder = Cylinder(center, axis, height, radius, interface)
        back_cap = SphericalCap(center, -axis, radius, pi / 2, interface)
        new{N}((front_cap, cylinder, back_cap))
    end
end
