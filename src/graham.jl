
function pseudo_angle(x :: T, y :: T) where {T <: Real}
    a = x // (abs(x) + abs(y))
    y < 0 && return a - 1
    return 1 - a
end
pseudo_angle(p :: Tuple{T,T}) where {T <: Real} = pseudo_angle(p[1], p[2])

ccw(a::Tuple{T,T}, b::Tuple{T,T}, c::Tuple{T,T}) where {T <: Real} =
(b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

function graham_scan!(points :: Vector{Tuple{T,T}}) where {T <: Real}
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

graham_scan(points :: Vector{Tuple{T,T}}) where {T <: Real} = graham_scan!(deepcopy(points))

convex_hull(points :: Vector{Tuple{T,T}}) where {T <: Real} =
RationalPolygon(graham_scan(points))
