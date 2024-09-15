module TraycerElements

using StaticArrays
using TraycerCore
import TraycerCore: origin, direction, normal, interface, tracingtypes, isobject, primitives
import TraycerCore: minintersection!, distance!, intensity, reset!, onintersect

function quadraticroots(a::Float64, b::Float64, c::Float64)::Tuple{Float64,Float64}
    temp = b^2 - 4 * a * c
    temp < 1e-10 && return 0.0, 0.0
    return (-b - sqrt(temp)) / 2a, (-b + sqrt(temp)) / 2a
end

include("surfaces/surfaces.jl")
include("objects/objects.jl")
include("compoundsurfaces/compoundsurfaces.jl")
include("compoundobjects/compoundobjects.jl")
include("coatedelements/coatedelement.jl")
include("linkedelements/linkedelement.jl")
end
