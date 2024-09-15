struct Bubble{N} <: CompoundObject{N}
    primitives::Tuple{Sphere{N},Sphere{N}}

    function Bubble(outer_sphere::Sphere{N}, inner_sphere::Sphere{N}) where {N}
        @assert radius(outer_sphere) > radius(inner_sphere) + 1e-10 "Outer sphere must have larger radius"
        new{N}((outer_sphere, inner_sphere))
    end
end

outer_sphere(bubble::Bubble) = bubble.primitives[1]
inner_sphere(bubble::Bubble) = bubble.primitives[2]
outerradius(bubble::Bubble) = radius(outer_sphere(bubble))
innerradius(bubble::Bubble) = radius(inner_sphere(bubble))
center(bubble::Bubble) = center(outer_sphere(bubble))

function Bubble(center::SVector{N,<:Real}, outerradius::Real, innerradius::Real, outer_refractive_index::Real, inner_refractive_index::Real=1.0) where {N}
    outer_sphere = Sphere(center, outerradius, outer_refractive_index)
    inner_sphere = Sphere(center, innerradius, inner_refractive_index)
    Bubble(outer_sphere, inner_sphere)
end

show(io::IO, bubble::Bubble) = print(io, "Bubble($(center(bubble)), $(outerradius(bubble)), $(innerradius(bubble)))")
