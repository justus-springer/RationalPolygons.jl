abstract type SliceBehaviour{T <: Integer} end

struct NoSlice{T <: Integer} <: SliceBehaviour{T} end

struct SliceThroughVertex{T <: Integer} <: SliceBehaviour{T}
    p :: RationalPoint{T}
end

struct SliceThroughBoundary{T <: Integer} <: SliceBehaviour{T}
    p :: RationalPoint{T}
    q :: RationalPoint{T}
end

struct SliceThroughInterior{T <: Integer} <: SliceBehaviour{T}
    p :: RationalPoint{T}
    q :: RationalPoint{T}
end

slice_length(::NoSlice) = 0
slice_length(::SliceThroughVertex) = 0
slice_length(sl :: SliceThroughBoundary) = distance(sl.p, sl.q)
slice_length(sl :: SliceThroughInterior) = distance(sl.p, sl.q)


@doc raw"""
    slice(L :: Line{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer}

Compute the intersection of the line `L` with the polygon `P`. Possible output
values are:

- `NoSlice()`, if the line does not intersect with the polygon,
- `SliceThroughVertex(p)`, if the line intersects `P` in precisely one vertex `p`,
- `SliceThroughBoundary(p,q)`, if the line intersects `P` in precisely two vertices `p` and `q` or
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
    return SliceThroughBoundary{T}(points[1],points[2])
end


@doc raw"""
    slice_length(L :: Line{T}, P :: RationalPolygon{T}) where {T <: Integer}
    
Return the length of the slice given by the intersection of `L` and `P`.
Returns zero if `L` and `P` do not intersect.

"""
slice_length(L :: Line{T}, P :: RationalPolygon{T}) where {T <: Integer} =
slice_length(slice(L,P))
