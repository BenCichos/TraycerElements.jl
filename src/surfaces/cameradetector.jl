mutable struct CameraDetector{N} <: Surface{N}
    const plane::Plane{N}
    const position::Function
    const counter::Matrix{Float64}
    const direction_fn::Function
    const angleresolution::Int
    state::Int

    function CameraDetector(origin::SVector{2,<:Real}, normal::SVector{2,<:Real}, size::Real; interface::Interface=ConsumingInterface(), resolution::Int=128, viewangle::Real=0.0, angleresolution::Int=1)
        plane = Plane(origin, normal, interface, size)
        rotation_matrix = rotation(normal)
        position = (i::Int) -> (rotation_matrix * SVector(i, 0.0)) + center(plane)
        counter = zeros(resolution, 1)
        new{2}(plane, position, counter, camera_ray_direction_fn(normal, viewangle), angleresolution)
    end

    function CameraDetector(origin::SVector{3,<:Real}, normal::SVector{3,<:Real}, size::Tuple{<:Real,<:Real}; interface::Interface=ConsumingInterface(), resolution::Int=128, viewangle::Real=0.0, angleresolution::Int=1)
        plane = Plane(origin, normal, interface, size)
        rotation_matrix = rotation(normal)
        position = (i::Int, j::Int) -> (rotation_matrix * (SVector((i - 0.5) * size[1] / resolution, (j - 0.5) * size[2] / resolution, 0.0) - SVector(size[1] / 2, size[2] / 2, 0.0))) + center(plane)
        counter = zeros(resolution, resolution)
        new{3}(plane, position, counter, camera_ray_direction_fn(normal, viewangle), angleresolution)
    end
end
export CameraDetector

plane(camera::CameraDetector) = camera.plane
center(camera::CameraDetector) = center(plane(camera))
normal(camera::CameraDetector) = normal(plane(camera))
size(camera::CameraDetector) = size(plane(camera))
length(camera::CameraDetector) = length(camera.counter)
state(camera::CameraDetector) = camera.state
direction_fn(camera::CameraDetector) = camera.direction_fn
angleresolution(camera::CameraDetector) = camera.angleresolution

resolution(camera::CameraDetector) = size(camera.counter)
position(camera::CameraDetector) = camera.position
counter(camera::CameraDetector) = camera.counter
reset!(camera::CameraDetector) = fill!(camera.counter, 0.0)
state!(camera::CameraDetector, new_state::Int) = camera.state = new_state

show(io::IO, camera::CameraDetector) = print(io, "CameraDetector($(center(camera)), $(normal(camera)), $(size(camera)), $(resolution(camera))")
BoundingBox(camera::CameraDetector) = BoundingBox(plane(camera))

function onintersect(camera::CameraDetector{3}, ray::Ray{3}, ::Float64, normal::SVector{3,<:Real})
    cos_with_normal = abs(dot(direction(ray), normal))
    cos_with_normal < 0.9 && return nothing
    counter(camera)[state(camera)] += intensity(ray) * cos_with_normal
    return nothing
end

function iterate(c::CameraDetector{3}, state::Int=1)
    state > length(c) && return nothing
    state!(c, state)
    i, j = fldmod1(state, resolution(c)[1])
    isone(angleresolution(c)) && return ((Ray(position(c)(i, j), normal(c)) for _ in 1:1), state + 1)
    interpolate_range = range(0, 1, angleresolution(c))
    gen = (
        Ray(position(c)(i, j), direction_fn(c)(c0, c1))
        for c0 in interpolate_range
        for c1 in interpolate_range
        if (iszero(c1) && iszero(c0) || !iszero(c1))
    )
    (gen, state + 1)
end

@inline perp(v::SVector{3,<:Real}) = SVector(v[3], -v[2], v[1])

camera_ray_direction_fn(normal::SVector{3,<:Real}, angle::Real) = ((c0::Real, c1::Real) -> rotate_vector(quatwithaxis(normal, 2pi * c0) * quatwithaxis(perp(normal), angle * c1), normal))
