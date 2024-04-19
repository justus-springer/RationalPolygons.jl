
struct RationalPolygon{T<:Integer}
    rationality :: T
    vs :: Vector{LatticePoint{T}}

    RationalPolygon(rationality :: T, vs :: Vector{LatticePoint{T}}) where {T <: Integer} =
    new{T}(rationality, vs)

    function RationalPolygon(points :: Vector{<:Point{T}}) where {T <: Integer}
        r = lcm(rationality.(points))
        integral_points = numerator.(r .* points)
        return new{T}(r,integral_points)
    end

end

Base.show(io :: IO, P :: RationalPolygon) =
Base.print(io, "Rational polygon of rationality $(P.rationality)")

Base.:(==)(P :: RationalPolygon, Q :: RationalPolygon) =
P.rationality == Q.rationality && P.vs == Q.vs

rationality(P :: RationalPolygon) = P.rationality

lattice_vertices(P :: RationalPolygon) = P.vs

vertices(P :: RationalPolygon) = lattice_vertices(P) .// rationality(P)
