
@doc raw"""
    abstract type RationalPolygon{T <: Integer} end

A rational polygon in two-dimensional rational space. Subtypes must at least
implement `vertices` and `affine_halfplanes`.

"""
abstract type RationalPolygon{T <: Integer} end

Base.in(x :: Point{T}, P :: RationalPolygon) where {T <: Integer} =
all(H -> x ∈ H, affine_halfplanes(P))

@attr number_of_vertices(P :: RationalPolygon) = length(vertices(P))

@attr rationality(P :: RationalPolygon) = lcm(rationality.(vertices(P)))

@attr function edges(P :: RationalPolygon)
    vs, r = vertices(P), number_of_vertices(P)
    return [(vs[i], vs[mod(i+1,1:r)]) for i = 1 : r]
end

Base.show(io :: IO, P :: RationalPolygon) =
Base.print(io, "Rational polygon of rationality $(rationality(P)) with $(number_of_vertices(P)) vertices.")
