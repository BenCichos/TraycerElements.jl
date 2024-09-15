include("lens.jl")
export Lens, front_cap, cylinder, back_cap, center, radius, thickness, radius_curvature, axis

include("bubble.jl")
export Bubble, outer_sphere, inner_sphere, outerradius, innerradius, center

include("closedcylinder.jl")
export ClosedCylinder, cylinder, front_plane, back_plane, center, radius, height

include("pill.jl")
export Pill, front_cap, cylinder, back_cap, center, radius
