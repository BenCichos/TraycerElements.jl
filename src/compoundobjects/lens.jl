struct Lens{N} <: CompoundObject{N}
    primitives::Tuple{SphericalCap{N},Cylinder{N},SphericalCap{N}}

    function Lens(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, cap_radius::Real, radius_lens::Real, thickness::Real, interface) where {N}
        @assert cap_radius >= radius_lens "The radius of the spherical caps must be larger than or equal to the radius of the the lens"

        axis = normalize(axis)
        theta = asin(radius_lens / cap_radius)
        cap_thickness = cap_radius*(1 - cos(theta))
        cap_offset = axis *(cap_radius - cap_thickness - thickness / 2)

        front_cap = SphericalCap(center - cap_offset, axis, cap_radius, theta, interface)
        back_cap = SphericalCap(center + cap_offset, -axis, cap_radius, theta, interface)
        cylinder = Cylinder(center - (thickness / 2 * axis), axis, thickness, radius_lens, interface)

        new{N}((front_cap, cylinder, back_cap))
    end
end

front_cap(lens::Lens) = lens.primitives[1]
cylinder(lens::Lens) = lens.primitives[2]
back_cap(lens::Lens) = lens.primitives[3]
center(lens::Lens) = center(cylinder(lens))
radius(lens::Lens) = radius(cylinder(lens))
thickness(lens::Lens) = height(cylinder(lens))
radius_curvature(lens::Lens) = radius(front_cap(lens))
axis(lens::Lens) = axis(cylinder(lens))

function Lens(center::SVector{N,<:Real}, axis::SVector{N,<:Real}, cap_radius::Real, radius_lens::Real, thickness::Real, refractive_index::Real) where {N}
    Lens(center, axis, cap_radius, radius_lens, thickness, Interface(refractive_index))
end

show(io::IO, lens::Lens) = print(io, "Lens($(center(lens)), $(axis(lens)), $(radius(lens)), $(radius_curvature(lens)), $(thickness(lens)))")
