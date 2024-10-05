struct ClosedCylinder{N} <: CompoundObject{N}
    primitives::Tuple{CircularPlane{N},Cylinder{N},CircularPlane{N}}

    function ClosedCylinder(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, radius::Real, height::Real, interface) where {N}
        front_plane = CircularPlane(center + height * axis, axis, radius, interface)
        cylinder = Cylinder(center, axis, height, radius, interface)
        back_plane = CircularPlane(center, -axis, radius, interface)
        new{N}((front_plane, cylinder, back_plane))
    end
end
