include("circularplane.jl")
export CircularPlane

include("plane.jl")
export Plane

include("sphericalcap.jl")
export SphericalCap

include("detector.jl")
export Detector

include("cylinder.jl")
export Cylinder

include("ring.jl")
export Ring

include("cameradetector.jl")
export CameraDetector, counter

include("lightplane.jl")
export LightPlane


include("lightcone.jl")
export LightCone

include("camera.jl")
export Camera, counter, reset!, state!
