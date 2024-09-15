struct Box{N} <: Object{N}
    limits::NTuple{N,Tuple{Float64,Float64}}
    interface::Interface

    function Box(limits::NTuple{2,Tuple{<:Real,<:Real}}, interface)
        new{2}(limits, interface)
    end

    function Box(limits::NTuple{3,Tuple{<:Real,<:Real}}, interface)
        new{3}(limits, interface)
    end
end

SystemBox(limits::NTuple{2,Tuple{<:Real,<:Real}}, m::Material) = Box(limits, ConsumingInterface(m))
SystemBox(limits::NTuple{3,Tuple{<:Real,<:Real}}, m::Material) = Box(limits, ConsumingInterface(m))
SystemBox(xlimits::Tuple{<:Real,<:Real}, ylimits::Tuple{<:Real,<:Real}, m::Material) = Box((xlimits, ylimits), m)
SystemBox(xlimits::Tuple{<:Real,<:Real}, ylimits::Tuple{<:Real,<:Real}, zlimits::Tuple{<:Real,<:Real}, m::Material) = SystemBox((xlimits, ylimits, zlimits), m)

limits(box::Box) = box.limits
xlimits(box::Box) = box.limits[1]
ylimits(box::Box) = box.limits[2]
zlimits(box::Box{3}) = box.limits[3]

show(io::IO, box::Box) = print(io, "Box($(limits(box)), $(interface(box)))")
function normal(box::Box{2}, point::SVector{2,Float64})
    point[1] in xlimits(box) && return SVector(-1.0 * sign(point[1]), 0.0)
    SVector(0.0, -1.0 * sign(point[2]))
end
function normal(box::Box{3}, point::SVector{3,Float64})
    (point[1] in xlimits(box)) && return SVector(-1.0 * sign(point[1]), 0.0, 0.0)
    (point[2] in ylimits(box)) && return SVector(0.0, -1.0 * sign(point[2]), 0.0)
    SVector(0.0, 0.0, -1.0 * sign(point[2]))
end

function minintersection!(minintersection::MinIntersection, box::Box{2}, ray::Ray{2})
    direction_x, direction_y = direction(ray)
    origin_x, origin_y = origin(ray)

    boundary_x = xlimits(box)[signbit(direction_x) ? 1 : 2]
    boundary_y = ylimits(box)[signbit(direction_y) ? 1 : 2]

    alpha_x = iszero(direction_x) ? Inf : (boundary_x - origin_x) / direction_x
    alpha_y = iszero(direction_y) ? Inf : (boundary_y - origin_y) / direction_y

    alpha = min(alpha_x, alpha_y)
    alpha <= 1e-10 && return nothing
    isinf(alpha) && return nothing
    distance!(minintersection, alpha)
end

function minintersection!(minintersection::MinIntersection, box::Box{3}, ray::Ray{3})
    direction_x, direction_y, direction_z = direction(ray)
    origin_x, origin_y, origin_z = origin(ray)

    boundary_x = xlimits(box)[signbit(direction_x) ? 1 : 2]
    boundary_y = ylimits(box)[signbit(direction_y) ? 1 : 2]
    boundary_z = zlimits(box)[signbit(direction_z) ? 1 : 2]

    alpha_x = iszero(direction_x) ? Inf : (boundary_x - origin_x) / direction_x
    alpha_y = iszero(direction_y) ? Inf : (boundary_y - origin_y) / direction_y
    alpha_z = iszero(direction_z) ? Inf : (boundary_z - origin_z) / direction_z

    alpha = min(alpha_x, alpha_y, alpha_z)
    alpha <= 1e-10 && return nothing
    isinf(alpha) && return nothing
    distance!(minintersection, alpha)
end
