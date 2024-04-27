
@doc raw"""
    abstract type RationalPolygon{T <: Integer} end

A rational polygon in two-dimensional rational space. Subtypes must at least
implement `rationality`, `vertices` and `affine_halfplanes`.

"""
abstract type RationalPolygon{T <: Integer} end

Base.in(x :: Point{T}, P :: RationalPolygon) where {T <: Integer} =
all(H -> x ∈ H, affine_halfplanes(P))

@attr minimal_rationality(P :: RationalPolygon) = lcm(rationality.(vertices(P)))

@attr number_of_vertices(P :: RationalPolygon) = length(vertices(P))

@attr lattice_vertices(P :: RationalPolygon) = numerator.(rationality(P) .* vertices(P))

@attr function edges(P :: RationalPolygon)
    vs, r = vertices(P), number_of_vertices(P)
    return [(vs[i], vs[mod(i+1,1:r)]) for i = 1 : r]
end

Base.show(io :: IO, P :: RationalPolygon) =
Base.print(io, "Rational polygon of rationality $(rationality(P)) with $(number_of_vertices(P)) vertices.")

@doc raw"""
    is_maximal(P :: RationalPolygon)

Check whether a rational polygon is maximal among all polygons sharing
the same rationality and number of interior lattice points.

"""
@attr function is_maximal(P :: RationalPolygon)
    k = rationality(P)
    ps = k_rational_points(k,P)
    Q = intersect_halfplanes(affine_halfplanes(P) .- 1 // k)
    for p ∈ boundary_k_rational_points(k, Q)
        new_P = convex_hull([ps ; p]; rationality = k)
        number_of_interior_lattice_points(new_P) <= number_of_interior_lattice_points(P) && return false
    end
    return true
end
