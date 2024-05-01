ccw(a::Tuple{T,T}, b::Tuple{T,T}, c::Tuple{T,T}) where {T <: Real} =
(b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

function graham_scan!(points :: Vector{Tuple{T,T}}) where {T <: Real}
    # remove duplicate points
    unique!(points)

    # handle edge cases
    length(points) < 3 && return points

    # find an extremal point, i.e. a point that is definetely a vertex
    # of the convex hull. We take the lexicographical minimum here.
    i0 = argmin(points)
    points[1], points[i0] = points[i0], points[1]
    p0 = points[1]

    # compute angle and distance with p0 for all remaining points
    #
    other_points_with_data = [(pseudo_angle_with_distance(p - p0), p) for p ∈ points[2:end]]
    # sort the remaining points lexicographically by their pair (angle, distance)
    sort!(other_points_with_data)

    # for collinear points, only keep the one furthest away from p0, since
    # the others are clearly not vertices. This is done by traversing
    # the sorted list from above and always adding the last point of all
    # sets of points sharing the same angle.
    empty!(points)
    push!(points, p0)
    for i = 1 : length(other_points_with_data)-1
        (a1, _), p1 = other_points_with_data[i]
        (a2, _), p2 = other_points_with_data[i+1]
        if a1 != a2
            push!(points, p1)
        end
    end
    push!(points, last(other_points_with_data)[2])

    # the main part of the algorithm
    n = length(points)
    k = 2
    for i = 3 : n
        while ccw(points[k-1], points[k], points[i]) <= 0
            if k > 2
                k -= 1
            else # (i == n)
                break
            end
        end
        k += 1
        points[i], points[k] = points[k], points[i]
    end

    return points[1:k]
end

graham_scan(points :: Vector{Tuple{T,T}}) where {T <: Real} = graham_scan!(deepcopy(points))
