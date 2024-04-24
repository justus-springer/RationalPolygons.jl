
@attributes mutable struct ConvexHull{T<:Integer} <: RationalPolygon{T}
    vertices :: Vector{RationalPoint{T}}

    ConvexHull(rationality :: T, vs :: Vector{LatticePoint{T}}) where {T <: Integer} =
    new{T}(rationality, vs)

    function ConvexHull(vertices :: Vector{<:Point{T}}) where {T <: Integer}
        return new{T}(vertices .// 1)
    end

end

Base.:(==)(P :: ConvexHull, Q :: ConvexHull) =
P.vertices == Q.vertices

vertices(P :: ConvexHull) = P.vertices

@attr affine_halfplanes(P :: ConvexHull) =
[AffineHalfplaneByLine(LineThroughPoints(e[1], e[2])) for e ∈ edges(P)]

convex_hull(points :: Vector{<:Point{T}}) where {T <: Integer} =
ConvexHull(graham_scan(points))

