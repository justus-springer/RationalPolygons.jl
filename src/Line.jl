
@doc raw"""
    abstract type Line{T <: Integer} end

A line in 2-dimensional rational space. Subtypes must at least implement:
`base_point` and `direction_vector`.

"""
abstract type Line{T <: Integer} end

function Base.in(x :: Point{T}, L :: Line{T}) where {T <: Integer}
    x0 = x - base_point(L)
    v = direction_vector(L)
    return iszero(det(v,x0))
end

point_by_parameter(L :: Line{T}, t :: T) where {T <: Integer} =
base_point(L) + t * direction_vector(L)

reverse_direction(L :: Line{T}) where {T <: Integer} =
LineByDirection(base_point(L), -direction_vector(L))

@doc raw"""
    struct LineByDirection{T<:Integer} <: Line{T}

A line described by a base point and direction vector.

"""
struct LineByDirection{T<:Integer} <: Line{T}
    A :: RationalPoint{T}
    dir :: RationalPoint{T}
    function LineByDirection(A :: Point{T}, dir :: Point{T}) where {T <: Integer}
        !iszero(dir) || error("direction vector can't be zero")
        new{T}(A,dir)
    end
end

direction_vector(L :: LineByDirection) = L.dir
base_point(L :: LineByDirection) = L.A

Base.show(io :: IO, L :: LineByDirection) =
print(io, "Line with base point $(L.A) and direction vector $(L.dir)")

@doc raw"""
    struct LineThroughPoints{T<:Integer} <: Line{T}

A line described by two points in rational 2D space.

"""
struct LineThroughPoints{T<:Integer} <: Line{T}
    A :: RationalPoint{T}
    B :: RationalPoint{T}
    function LineThroughPoints(A :: Point{T}, B :: Point{T}) where {T <: Integer}
        A ≠ B || error("points cannot be equal")
        new{T}(A,B)
    end
end

direction_vector(L :: LineThroughPoints) = L.B - L.A
base_point(L :: LineThroughPoints) = L.A

Base.show(io :: IO, L :: LineThroughPoints) =
print(io, "Line through the points $(L.A) and $(L.B)")

struct HorizontalLine{T <: Integer} <: Line{T}
    y :: Rational{T}
    HorizontalLine(y :: Rational{T}) where {T <: Integer} = new{T}(y)
    HorizontalLine(y :: T) where {T <: Integer} = new{T}(y // 1)
end

base_point(L :: HorizontalLine) = (zero(L.y), L.y)
direction_vector(L :: HorizontalLine) = (one(L.y), zero(L.y))

Base.show(io :: IO, L :: HorizontalLine) =
print(io, "Horizontal line at y = $(L.y)")

struct VerticalLine{T <: Integer} <: Line{T}
    x :: Rational{T}
    VerticalLine(y :: Rational{T}) where {T <: Integer} = new{T}(y)
    VerticalLine(y :: T) where {T <: Integer} = new{T}(y // 1)
end

base_point(L :: VerticalLine) = (L.x, zero(L.x))
direction_vector(L :: VerticalLine) = (zero(L.x), one(L.x))

Base.show(io :: IO, L :: VerticalLine) =
print(io, "Vertical line at y = $(L.x)")


abstract type IntersectionBehaviour end

struct IntersectInPoint{T <: Integer} <: IntersectionBehaviour
    p :: Point{T}
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

    if det(v1,v2) == 0
        return x1 ∈ L2 ? LinesAreEqual() : NoIntersection()
    end

    A = QQ[v1[1] -v2[1] ; v1[2] -v2[2]]
    b = QQ[x2[1]-x1[1] ; x2[2]-x1[2]]

    t = Rational{T}(can_solve_with_solution(A,b; side = :right)[2][1])
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
    k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}) where {T <: Integer}
    
Return all `k`-rational points lying on the line segment from p to q,
*excluding* p and q themselves.

"""
function k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}) where {T <: Integer}

    if q[1] < p[1] || (p[1] == q[1] && q[2] < p[2])
        return integral_primitive_points_on_line_segment(k,q,p)
    end

    L = LineThroughPoints(p,q)
    res = RationalPoint{T}[]
    if p[1] < q[1]
        for t = floor_k_rational(k, p[1] + 1 // k) : 1 // k : ceil_k_rational(k, q[1] - 1 // k)
            x = intersection_point(L, VerticalLine(t))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    elseif p[2] < q[2]
        for t = floor_k_rational(k, p[2] + 1 // k) : 1 // k : ceil_k_rational(k, q[2] - 1 // k)
            x = intersection_point(L, HorizontalLine(t))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    end

    return res

end

integral_points_on_line_segment(p :: Point{T}, q :: Point{T}) where {T <: Integer} =
k_rational_points_on_line_segment(1,p,q)

@doc raw"""
    next_k_rational_point(p :: Point{T}, L :: Line{T}) where {T <: Integer}

Given a point `p` and a line `L` with `p ∈ L`, return the next `k`-rational
point on `L` after `p`, following the direction vector of the line, excluding
`p` itself.

"""
function next_k_rational_point(k :: T, p :: Point{T}, L :: Line{T}) where {T <: Integer}
    p ∈ L || error("point must lie on the line")
    v = direction_vector(L)
    if v[1] ≠ 0
        t = v[1] > 0 ? floor_k_rational(k, p[1]+1//k) : ceil_k_rational(k,p[1]-1//k)
        x = intersection_point(L, VerticalLine(t))
        while !is_k_rational(k,x)
            t += sign(v[1]) // k
            x = intersection_point(L, VerticalLine(t))
        end
        return x
    elseif v[2] ≠ 0
        t = v[2] > 0 ? floor_k_rational(k, p[2]+1//k) : ceil_k_rational(k, p[2]-1//k)
        x = intersection_point(L, HorizontalLine(t))
        while !is_k_rational(k,x)
            t += sign(v[2]) // k
            x = intersection_point(L, HorizontalLine(t))
        end
        return x
    end
end

next_integral_point(p :: Point{T}, L :: Line{T}) where {T <: Integer} =
next_k_rational_point(1, p, L)

@doc raw"""
    previous_k_rational_point(k :: T, p :: Point{T}, L :: Line{T}) where {T <: Integer}

Given a point `p` and a line `L` with `p ∈ L`, return the previous integral
primitive point on `L` after `p`, following the direction vector of the line,
excluding `p` itself.

"""
previous_k_rational_point(k :: T, p :: Point{T}, L :: Line{T}) where {T <: Integer} =
next_k_rational_point(k, p, reverse_direction(L))

previous_integral_point(p :: Point{T}, L :: Line{T}) where {T <: Integer} =
previous_k_rational_point(1, p, L)
