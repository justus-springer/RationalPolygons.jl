@doc raw"""
    SliceBehaviour{T <: Integer}

An abstract supertype for possible intersection behaviours of a line with a
polygon. There are four subtypes `NoSlice`, `SliceThroughVertex`,
`SliceThroughEdge` and `SliceThroughInterior`.

"""
abstract type SliceBehaviour{T <: Integer} end


@doc raw"""
    NoSlice{T <: Integer}

The intersection behaviour of a line that does not intersect with a polygon.
This struct has no fields.

"""
struct NoSlice{T <: Integer} <: SliceBehaviour{T} end


@doc raw"""
    SliceThroughVertex{T <: Integer}

The intersection behaviour of a line that intersects a polygon in a single
vertex. This struct has a single field `p`, which is the vertex.

"""
struct SliceThroughVertex{T <: Integer} <: SliceBehaviour{T}
    p :: RationalPoint{T}
end


@doc raw"""
    SliceThroughEdge{T <: Integer}

The intersection behaviour of a line that intersects a polygon in two adjacent
vertices. This struct has two fields `p` and `q`, which are the vertices that
the line intersects with.

"""
struct SliceThroughEdge{T <: Integer} <: SliceBehaviour{T}
    p :: RationalPoint{T}
    q :: RationalPoint{T}
end


@doc raw"""
    SliceThroughInterior{T <: Integer}

The intersection behaviour of a line that intersects a polyon non-trivially in
its interior. This struct has two fields `p` and `q`, which are the
intersection points of the line with the boundary of the polygon.

"""
struct SliceThroughInterior{T <: Integer} <: SliceBehaviour{T}
    p :: RationalPoint{T}
    q :: RationalPoint{T}
end


slice_length(::NoSlice) = 0
slice_length(::SliceThroughVertex) = 0
slice_length(sl :: SliceThroughEdge) = distance(sl.p, sl.q)
slice_length(sl :: SliceThroughInterior) = distance(sl.p, sl.q)


@doc raw"""
    slice(L :: Line{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer}

Compute the intersection of the line `L` with the polygon `P`. Possible output
values are:

- `NoSlice()`, if the line does not intersect with the polygon,
- `SliceThroughVertex(p)`, if the line intersects `P` in precisely one vertex `p`,
- `SliceThroughEdge(p,q)`, if the line intersects `P` in two adjacent vertices `p` and `q` or
- `SliceThroughInterior(p,q)`, if the line passes through the interior of `P` and intersects the boundary of `P` in precisely the points `p` and `q`.

"""
function slice(L :: Line{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer}
    points = RationalPoint{T}[]
    through_interior = true
    for i = 1 : N
        int_behaviour = intersection_behaviour(L, line_through_points(P[i], P[i+1]))
        if int_behaviour isa IntersectInPoint{T}
            p = int_behaviour.p
            p ∈ P || continue
            push!(points, p)
        elseif int_behaviour isa LinesAreEqual{T}
            through_interior = false
        end
    end

    isempty(points) && return NoSlice{T}()
    unique!(points)
    length(points) == 1 && return SliceThroughVertex{T}(points[1])
    through_interior && return SliceThroughInterior{T}(points[1],points[2])
    return SliceThroughEdge{T}(points[1],points[2])
end


@doc raw"""
    slice_length(L :: Line{T}, P :: RationalPolygon{T}) where {T <: Integer}
    
Return the length of the slice given by the intersection of `L` and `P`.
Returns zero if `L` and `P` do not intersect.

"""
slice_length(L :: Line{T}, P :: RationalPolygon{T}) where {T <: Integer} =
slice_length(slice(L,P))
