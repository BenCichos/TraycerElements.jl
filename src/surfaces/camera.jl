mutable struct Camera{N} <: Surface{N}
    const plane::Plane{N}
    const position::Function
    const counter::Matrix{Float64}
    const wavelength_range::LinRange{Float64,Int64}
    const ray_index::Float64
    const iterator_length::Int
    state::Int

    function Camera(origin::SVector{3,<:Real}, normal::SVector{3,<:Real}, focal_length::Real, size::Tuple{<:Real,<:Real}; wavelength_range::Tuple{Float64,Float64}=(0.54, 0.54), wavelength_resolution::Int=1, ray_index::Float64=1.0, interface::Interface=ConsumingInterface(), resolution::Int=128)
        new_plane_origin = origin - focal_length * normal
        plane = Plane(new_plane_origin, normal, interface, size)
        rotation_matrix = rotation(normal)
        position = (i::Int, j::Int) -> (rotation_matrix * (SVector((i - 0.5) * size[1] / resolution, (j - 0.5) * size[2] / resolution, 0.0) - SVector(0.5 * size[1], 0.5 * size[2], -focal_length))) + center(plane)
        counter = zeros(resolution, resolution)
        new{3}(plane, position, counter, LinRange(wavelength_range..., wavelength_resolution), ray_index, length(counter) * wavelength_resolution, 1)
    end
end
export Camera

plane(camera::Camera) = camera.plane
center(camera::Camera) = center(plane(camera))
normal(camera::Camera) = normal(plane(camera))
size(camera::Camera) = size(plane(camera))
counter_length(camera::Camera) = length(counter(camera))
length(camera::Camera) = camera.iterator_length
ray_index(camera::Camera) = camera.ray_index
state(camera::Camera) = camera.state

resolution(camera::Camera) = size(camera.counter)
position(camera::Camera) = camera.position
counter(camera::Camera) = camera.counter
wavelength_range(camera::Camera) = camera.wavelength_range
reset!(camera::Camera) = fill!(camera.counter, 0.0)
state!(camera::Camera, new_state::Int) = camera.state = new_state

show(io::IO, camera::Camera) = print(io, "Camera($(center(camera)), $(normal(camera)), $(size(camera)), $(resolution(camera))")
BoundingBox(camera::Camera) = BoundingBox(plane(camera))

function onintersect(camera::Camera{3}, ray::Ray{3}, ::Float64, normal::SVector{3,<:Real})
    counter(camera)[mod1(state(camera), counter_length(camera))] += intensity(ray) * abs(dot(direction(ray), normal))
    return nothing
end

function iterate(c::Camera{3}, state::Int=1)
    state > length(c) && return nothing
    state!(c, state)
    i, j = fldmod1(mod1(state, counter_length(c)), resolution(c)[2])
    (Ray(position(c)(i, j), normalize(position(c)(i, j) - center(c)), 1.0, wavelength_range(c)[fld1(state, counter_length(c))], ray_index(c)), state + 1)
end

function cameratrace(os::OpticalSystem{N}, camera::Camera{N}; depth::Int=10, threshold_intensity::Float64=1e-5, keeptracedrays::Bool=false) where {N}
    @assert depth >= 0 "depth must be a non-negative integer"
    @assert threshold_intensity >= 0 "threshold_intensity must be a non-negative number"

    recursive_tracing_vector = RecursiveTracingVector(2^3 * length(rays(os)), AbstractOpticalElement{N}[tracingtypes(os)...]; keeptracedrays=keeptracedrays)
    minintersection = MinIntersection()

    for pixel_ray in camera
        _recursivetrace!(recursive_tracing_vector, minintersection, pixel_ray, depth, threshold_intensity)
    end

    return nothing
end
export cameratrace
