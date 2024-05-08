
@doc raw"""
    Line{T <: Integer}

A line in 2-dimensional rational space.

"""
struct Line{T <: Integer}
    base_point :: RationalPoint{T}
    direction_vector :: RationalPoint{T}
    function Line(base_point :: Point{T}, direction_vector :: Point{T}) where {T <: Integer}
        !iszero(direction_vector) || error("direction vector can't be zero")
        new{T}(base_point, direction_vector)
    end
end

base_point(L :: Line{T}) where {T <: Integer} = L.base_point
direction_vector(L :: Line{T}) where {T <: Integer} = L.direction_vector

Base.show(io :: IO, L :: Line) =
print(io, "Line with base point $(base_point(L)) and direction vector $(direction_vector(L))")

Base.:(==)(L1 :: Line, L2 :: Line) =
base_point(L1) == base_point(L2) && direction_vector(L1) == direction_vector(L2)

function Base.hash(L :: Line, h :: UInt64) 
    h = hash(L.base_point, h)
    h = hash(L.direction_vector, h)
    return h
end

Base.in(x :: Point{T}, L :: Line{T}) where {T <: Integer} =
iszero(det(x - base_point(L), direction_vector(L)))

point_by_parameter(L :: Line{T}, t :: T) where {T <: Integer} =
base_point(L) + t * direction_vector(L)

reverse_direction(L :: Line{T}) where {T <: Integer} =
Line(base_point(L), -direction_vector(L))

function normal_vector(L :: Line{T}) where {T <: Integer}
    v = primitivize(direction_vector(L))
    return LatticePoint(-v[2], v[1])
end

line_through_points(A :: Point{T}, B :: Point{T}) where {T <: Integer} =
Line(A, B - A)

horizontal_line(y :: Union{T, Rational{T}}) where {T <: Integer} =
Line(Point(zero(y), y), Point(one(y), zero(y)))

vertical_line(x :: Union{T, Rational{T}}) where {T <: Integer} =
Line(Point(x, zero(x)), Point(zero(x), one(x)))



abstract type IntersectionBehaviour end

struct IntersectInPoint{T <: Integer} <: IntersectionBehaviour
    p :: RationalPoint{T}
end

struct NoIntersection <: IntersectionBehaviour end

struct LinesAreEqual <: IntersectionBehaviour end 


@doc raw"""
    intersection_behaviour(L1 :: Line{T}, L2 :: Line{T}) where {T <: Integer}

Given two lines in 2D rational space, return the intersection behaviour of the
two lines: Possible values are `LinesAreEqual()`, `NoIntersection()` and
`IntersectInPoint(p)` where `p` is the unique intersection point.


"""
function intersection_behaviour(L1 :: Line{T}, L2 :: Line{T}) where {T <: Integer}
    v1, v2 = direction_vector(L1), direction_vector(L2)
    x1, x2 = base_point(L1), base_point(L2)

    d = det(v1,-v2)
    if d == 0
        return x1 ∈ L2 ? LinesAreEqual() : NoIntersection()
    end

    t = (v2[1]*(x2[2]-x1[2]) - v2[2]*(x2[1]-x1[1])) // d
    p = x1 + t * v1

    return IntersectInPoint{T}(p)
end

@doc raw"""
    intersection_point(L1 :: Line{T}, L2 :: Line{T}) where {T <: RationalUnion}

Return the intersection point of two lines in 2D rational space. Throws an
error if the lines do not intersect uniquely.

"""
function intersection_point(L1 :: Line{T}, L2 :: Line{T}) where {T <: Integer}
    int = intersection_behaviour(L1,L2)
    int isa IntersectInPoint || error("the lines do not intersect uniquely")
    return int.p
end


@doc raw"""
    k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}
    
Return all `k`-rational points lying on the line segment from p to q.
If `interior` is set to false, this includes `p` and `q` themselves, if
they are `k`-rational. Otherwise, `p` and `q` are always excluded.

"""
function k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}

    if q[1] < p[1] || (p[1] == q[1] && q[2] < p[2])
        return k_rational_points_on_line_segment(k, q, p; interior)
    end

    L = line_through_points(p,q)
    res = RationalPoint{T}[]
    if p[1] < q[1]
        lower_bound = interior ? floor_k_rational(k, p[1]+1//k) : ceil_k_rational(k, p[1]) 
        upper_bound = interior ? ceil_k_rational(k, q[1]-1//k) : floor_k_rational(k, q[1])
        for t = lower_bound : 1 // k : upper_bound
            x = intersection_point(L, vertical_line(t))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    elseif p[2] < q[2]
        lower_bound = interior ? floor_k_rational(k, p[2]+1//k) : ceil_k_rational(k, p[2])
        upper_bound = interior ? ceil_k_rational(k, q[2]-1//k) : floor_k_rational(k, q[2])
        for t = lower_bound : 1 // k : upper_bound
            x = intersection_point(L, horizontal_line(t))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    end

    return res

end

integral_points_on_line_segment(p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer} =
k_rational_points_on_line_segment(1, p, q; interior)
