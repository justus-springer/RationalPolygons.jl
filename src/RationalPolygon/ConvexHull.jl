
@attributes mutable struct ConvexHull{T<:Integer} <: RationalPolygon{T}
    rationality :: T
    vertices :: Vector{RationalPoint{T}}

    function ConvexHull(vertices :: Vector{<:Point{T}}; rationality :: Union{Missing,T} = missing) where {T <: Integer}
        if ismissing(rationality)
            k = lcm(RationalPolygons.rationality.(vertices))
        else
            k = rationality
        end
        return new{T}(k, vertices)
    end

end

Base.:(==)(P :: ConvexHull, Q :: ConvexHull) =
P.rationality == Q.rationality && P.vertices == Q.vertices

rationality(P :: ConvexHull) = P.rationality

vertices(P :: ConvexHull) = P.vertices

@attr affine_halfplanes(P :: ConvexHull) =
[AffineHalfplaneByLine(LineThroughPoints(e[1], e[2])) for e ∈ edges(P)]

convex_hull(points :: Vector{<:Point{T}}; rationality :: Union{Missing,T} = missing) where {T <: Integer} =
ConvexHull(graham_scan(points); rationality)

