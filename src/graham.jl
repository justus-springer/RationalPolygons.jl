
function pseudo_angle(p :: Point{T}) where {T <: Integer}
    x,y = p[1],p[2]
    (x == 0 && y == 0) && return 0
    a = x // (abs(x) + abs(y))
    y < 0 && return a - 1
    return 1 - a
end

ccw(a::Point{T}, b::Point{T}, c::Point{T}) where {T <: Integer} =
(b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

function graham_scan!(points :: Vector{<:Point{T}}) where {T <: Integer}
    # remove duplicate points
    unique!(points)

    # find an extremal point, i.e. a point that is definetely a vertex
    # of the convex hull. We take the lexicographical minimum here.
    i0 = argmin(points)
    points[1], points[i0] = points[i0], points[1]
    p0 = points[1]

    # find elements sharing the same angle with p0 and only keep
    # the furthest one, since the others are clearly not vertices.
    redundant_points = []
    for p ∈ points[2:end]
        p ∈ redundant_points && continue
        collinear_points = filter(q -> pseudo_angle(p - p0) == pseudo_angle(q - p0), points[2:end])
        sort!(collinear_points, by = q -> distance(p0,q))
        append!(redundant_points, collinear_points[1:end-1])
    end
    filter!(p -> p ∉ redundant_points, points)

    # sort the remaining points by their angle with p0
    points[2:end] = sort(points[2:end], by = p -> pseudo_angle(p - p0))

    # the main part of the algorithm
    n = length(points)
    k = 2
    for i = 1 : n
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
