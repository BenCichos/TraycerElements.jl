module TraycerElements

using StaticArrays: SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using TraycerCore
using Quaternions

using Base: @kwdef
import Base: size

const _TOLERANCE = 1e-10

include("utils.jl")
include("surfaces/surfaces.jl")
include("objects/objects.jl")
include("compoundsurfaces/compoundsurfaces.jl")
include("compoundobjects/compoundobjects.jl")
include("coatedelements/coatedelement.jl")
include("linkedelements/linkedelement.jl")
end
