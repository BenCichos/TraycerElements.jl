struct Detector{N} <: Surface{N}
    plane::Plane{N}
    binsize::Tuple{Float64,Float64}
    counter::Matrix{Float64}
    isdirectional::Bool

    function Detector(origin::SVector{2,<:Real}, normal::SVector{2,<:Real}, size::Real, interface::Interface=NULL_INTERFACE; isdirectional::Bool=true, resolution::Int=128)
        plane = Plane(origin, normal, interface, size)
        binsize = (size / resolution, 0.0)
        counter = zeros(resolution, 1)
        new{2}(plane, binsize, counter, isdirectional)
    end

    function Detector(origin::SVector{3,<:Real}, normal::SVector{3,<:Real}, size::Tuple{<:Real,<:Real}, interface::Interface=NULL_INTERFACE; isdirectional::Bool=true, resolution::Int=128)
        plane = Plane(origin, normal, interface, size)
        binsize = size ./ resolution
        counter = zeros(resolution, resolution)
        new{3}(plane, binsize, counter, isdirectional)
    end
end

plane(detector::Detector) = detector.plane
center(detector::Detector) = center(plane(detector))
normal(detector::Detector) = normal(plane(detector))
size(detector::Detector{3}) = size(plane(detector))
isdirectional(detector::Detector) = detector.isdirectional

resolution(detector::Detector) = size(detector.counter)
binsize(detector::Detector) = detector.binsize
counter(detector::Detector) = detector.counter
reset!(detector::Detector) = fill!(detector.counter, 0.0)

function onintersect(detector::Detector{2}, ray::Ray{2}, distance::Float64, normal::SVector{2,Float64})
    detector_size = only(size(detector))
    detector_binsize = only(binsize(detector))
    detector_vector = SVector(normal[2], -normal[1])
    position_on_detector = dot(origin(ray, distance) - center(detector), detector_vector)
    limit = detector_size / 2
    abs(position_on_detector) < limit || return nothing
    binindex = div(position_on_detector + limit, detector_binsize) |> Int
    @inbounds counter(detector)[binindex+1] += intensity(ray)
    return nothing
end

function onintersect(detector::Detector{3}, ray::Ray{3}, distance::Float64, normal::SVector{3,Float64})
    x_size, y_size = size(detector)
    x_binsize, y_binsize = binsize(detector)
    rotation_matrix = rotation(normal)
    position_on_detector = inv(rotation_matrix) * (origin(ray, distance) - center(detector))
    x_limit, y_limit = x_size / 2, y_size / 2
    abs(position_on_detector[1]) < x_limit || return nothing
    abs(position_on_detector[2]) < y_limit || return nothing
    x_binindex = div(position_on_detector[1] + x_limit, x_binsize) |> Int
    y_binindex = div(position_on_detector[2] + y_limit, y_binsize) |> Int
    @inbounds counter(detector)[x_binindex+1, y_binindex+1] += intensity(ray)
    return nothing
end

show(io::IO, detector::Detector) = print(io, "Detector($(center(detector)), $(normal(detector)), $(size(detector)), $(resolution(detector))")


minintersection!(minintersection::MinIntersection, detector::Detector{N}, ray::Ray{N}) where {N} = minintersection!(minintersection, plane(detector), ray)
