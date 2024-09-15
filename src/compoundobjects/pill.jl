struct Pill{N} <: CompoundObject{N}
    primitives::Tuple{SphericalCap{N},Cylinder{N},SphericalCap{N}}

    function Pill(front_cap::SphericalCap{N}, cylinder::Cylinder{N}, back_cap::CircularPlane{N}) where {N}
        @assert dot(normal(front_cap), normal(back_cap)) == -1.0 "The normal vectors of the front and back planes must be antiparallel"
        @assert dot(normal(front_cap), axis(cylinder)) == 1.0 "The normal vector of the front plane must be parallel to the cylinder axis"
        @assert radius(front_cap) == radius(cylinder) == radius(back_cap) "The radii of the front plane, cylinder, and back plane must be equal"
        @assert center(front_cap) == center(cylinder) + height(cylinder) * axis(cylinder) "The center of the front plane must be at the top of the cylinder"
        @assert center(back_cap) == center(cylinder) "The center of the back plane must be at the bottom of the cylinder"
        new{N}((front_cap, cylinder, back_cap))
    end
end

front_cap(pill::Pill) = pill.primitives[1]
cylinder(pill::Pill) = pill.primitives[2]
back_cap(pill::Pill) = pill.primitives[3]
center(pill::Pill) = center(cylinder(pill))
radius(pill::Pill) = radius(cylinder(pill))


function Pill(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, radius::Real, height::Real, refractive_index::Real) where {N}
    front_cap = HalfSphere(center + height * axis, axis, radius, Interface(refractive_index))
    cylinder = Cylinder(center, axis, height, radius, Interface(refractive_index))
    back_cap = HalfSphere(center, -axis, radius, Interface(refractive_index))
    Pill(front_cap, cylinder, back_cap)
end

show(io::IO, pill::Pill) = print(io, "Pill($(center(pill)), $(axis(pill)), $(radius(pill)), $(height(pill)))")
