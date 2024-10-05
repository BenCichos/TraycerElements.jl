module TraycerElements

using TraycerCore
using TraycerCore: RecursiveTracingVector, _recursivetrace!
using TraycerCore: @approx
using StaticArrays: SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using Quaternions

import TraycerCore: primitives

using Base: @kwdef
import Base: size

include("surfaces/surfaces.jl")
include("objects/objects.jl")
include("compoundsurfaces/compoundsurfaces.jl")
include("compoundobjects/compoundobjects.jl")
include("coatedelements/coatedelement.jl")
include("linkedelements/linkedelement.jl")
end
