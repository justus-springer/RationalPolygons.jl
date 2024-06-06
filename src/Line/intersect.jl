
@doc raw"""
    IntersectionBehaviour{T <: Integer}

An abstract supertype for possible intersection behaviours of two lines.
There are the three subtypes `IntersectInPoint`, `NoIntersection` and
`LinesAreEqual`.

"""
abstract type IntersectionBehaviour{T <: Integer} end


@doc raw"""
    IntersectInPoint{T <: Integer}

The intersection behaviour of two lines intersecting in a unique point. This
struct has a single field `p`, which is the intersection point.

"""
struct IntersectInPoint{T <: Integer} <: IntersectionBehaviour{T}
    p :: RationalPoint{T}
end


@doc raw"""
    NoIntersection{T <: Integer}
    
The intersection behaviour of parallel lines.

"""
struct NoIntersection{T <: Integer} <: IntersectionBehaviour{T} end


@doc raw"""
    LinesAreEqual{T <: Integer}

The intersection behaviour of two equal lines.

"""
struct LinesAreEqual{T <: Integer} <: IntersectionBehaviour{T} end 


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
        return x1 ∈ L2 ? LinesAreEqual{T}() : NoIntersection{T}()
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
