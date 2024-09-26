module TraycerElements

using TraycerCore
using TraycerCore: RecursiveTracingVector, _recursivetrace!
using ApproximateRelations
using StaticArrays: SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using Quaternions

using Base: @kwdef
import Base: size

const _TOLERANCE = 1e-10

include("surfaces/surfaces.jl")
include("objects/objects.jl")
include("compoundsurfaces/compoundsurfaces.jl")
include("compoundobjects/compoundobjects.jl")
include("coatedelements/coatedelement.jl")
include("linkedelements/linkedelement.jl")
end
