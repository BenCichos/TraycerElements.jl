struct ClosedCylinder{N} <: CompoundObject{N}
    primitives::Tuple{CircularPlane{N},Cylinder{N},CircularPlane{N}}

    function ClosedCylinder(front_plane::CircularPlane{N}, cylinder::Cylinder{N}, back_plane::CircularPlane{N}) where {N}
        @assert dot(normal(front_plane), normal(back_plane)) == -1.0 "The normal vectors of the front and back planes must be antiparallel"
        @assert dot(normal(front_plane), axis(cylinder)) == 1.0 "The normal vector of the front plane must be parallel to the cylinder axis"
        @assert radius(front_plane) == radius(cylinder) == radius(back_plane) "The radii of the front plane, cylinder, and back plane must be equal"
        @assert center(front_plane) == center(cylinder) + height(cylinder) * axis(cylinder) "The center of the front plane must be at the top of the cylinder"
        @assert center(back_plane) == center(cylinder) "The center of the back plane must be at the bottom of the cylinder"
        new{N}((front_plane, cylinder, back_plane))
    end
end

front_plane(closed_cylinder::ClosedCylinder) = closed_cylinder.primitives[1]
cylinder(closed_cylinder::ClosedCylinder) = closed_cylinder.primitives[2]
back_plane(closed_cylinder::ClosedCylinder) = closed_cylinder.primitives[3]
axis(closed_cylinder::ClosedCylinder) = axis(cylinder(closed_cylinder))
center(closed_cylinder::ClosedCylinder) = center(cylinder(closed_cylinder))
radius(closed_cylinder::ClosedCylinder) = radius(cylinder(closed_cylinder))
height(closed_cylinder::ClosedCylinder) = height(cylinder(closed_cylinder))


function ClosedCylinder(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, radius::Real, height::Real, refractive_index::Real) where {N}
    front_plane = CircularPlane(center + height * axis, axis, radius, Interface(refractive_index))
    cylinder = Cylinder(center, axis, height, radius, Interface(refractive_index))
    back_plane = CircularPlane(center, -axis, radius, Interface(refractive_index))
    ClosedCylinder(front_plane, cylinder, back_plane)
end

show(io::IO, closed_cylinder::ClosedCylinder) = print(io, "ClosedCylinder($(center(closed_cylinder)), $(axis(closed_cylinder)), $(radius(closed_cylinder)), $(height(closed_cylinder)))")
