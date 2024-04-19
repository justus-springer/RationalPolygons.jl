
function pseudo_angle(p :: Point{T}) where {T <: Integer}
    x,y = p[1],p[2]
    a = x // (abs(x) + abs(y))
    y < 0 && return a - 1
    return 1 - a
end

ccw(a::Point{T}, b::Point{T}, c::Point{T}) where {T <: Integer} =
(b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

function graham_scan!(points :: Vector{<:Point{T}}) where {T <: Integer}
    n = length(points)

    i = argmin(points)
    points[i], points[1] = points[1], points[i]
    p0 = points[1]

    points[2:end] = sort(points[2:end], by = p -> pseudo_angle(p0 .- p))

    k = 2
    for i = 1:n
        while (ccw(points[k-1], points[k], points[i]) <= 0)
            if k > 2
                k -= 1
            elseif (i == n)
                break
            else
                i += 1
            end
        end

        k += 1
        points[i], points[k] = points[k], points[i]
    end

    return points[1:k]
end

graham_scan(points :: Vector{<:Point{T}}) where {T <: Integer} = graham_scan!(deepcopy(points))

convex_hull(points :: Vector{<:Point{T}}) where {T <: Integer} =
RationalPolygon(graham_scan(points))
