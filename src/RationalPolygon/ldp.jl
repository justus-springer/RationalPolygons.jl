@doc raw"""
    contains_origin_in_interior(P :: RationalPolygon)

Check whether `P` contains the origin in its interior.

"""
contains_origin_in_interior(P :: RationalPolygon{T}) where {T <: Integer} =
contains_in_interior(LatticePoint{T}(0,0), P)


@doc raw"""
    is_primitive(P :: RationalPolygon)

Check whether all vertices of the scaled lattice polygon `rationality(P) * P`
are primitive.

"""
is_primitive(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(i -> is_primitive(scaled_vertex(P,i)), 1 : N)


@doc raw"""
    is_ldp(P :: RationalPolygon)

Check whether `P` is a ldp polygon, i.e. is primitive and contains the origin
in its interior.

"""
is_ldp(P :: RationalPolygon) = contains_origin_in_interior(P) && is_primitive(P)


function is_special_facet(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
    H1 = affine_halfplane(LatticePoint{T}(0,0), scaled_vertex(P,i))
    H2 = affine_halfplane(scaled_vertex(P,i+1), LatticePoint{T}(0,0))
    v = sum(vertices(P))
    return v ∈ H1 && v ∈ H2
end

special_facets(P :: RationalPolygon{T,N}) where {N, T <: Integer} =
filter(i -> is_special_facet(P,i), 1 : N)


@doc raw"""
    dual(P :: RationalPolygon{T}) where {T <: Integer}

Return the dual of a polygon `P`. Note that `P` must contain the origin in its interior.

"""
function dual(P :: RationalPolygon{T}) where {T <: Integer}
    contains_origin_in_interior(P) || error("this polygon does not contain the origin in its interior")
    Hs = affine_halfplanes(P)
    return convex_hull([normal_vector(H) // translation(H) for H ∈ Hs])
end


@doc raw"""
    multiplicity(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

The index of the sublattice spanned by the `i`-th and `i+1`-th scaled vertex of
`P` (i.e. the determinant of those two vertices). For ldp polygons, this equals
the order of the local class group associated with the toric fixed point
associated to the `i`-th and `i+1`-th ray.

"""
multiplicity(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} =
det(scaled_vertex(P,i), scaled_vertex(P,i+1))


@doc raw"""
    multiplicity(P :: RationalPolygon{T,N}) where {N,T <: Integer}

The order of the sublattice spanned by the scaled vertices of `P`. For ldp
polygons, this equals the order of the torsion part of the divisors class group
of the associated toric surface.

"""
multiplicity(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
gcd([multiplicity(P,i) for i = 1 : N])


@doc raw"""
    is_smooth(P :: RationalPolygon, i :: Int)

Check whether the cone spanned by the `i`-th and `i+1`-th vertex of `P` is
regular, i.e. generates the entire lattice. For ldp polygons, this means that
the toric fixed point associated to the `i`-th and `i+1`-th ray is smooth.

"""
is_smooth(P :: RationalPolygon, i :: Int) =
multiplicity(P, i) == 1


@doc raw"""
    is_smooth(P :: RationalPolygon)

Check whether all cones of the face fan of `P` are regular. For ldp polygons,
this means that the associated toric surface is smooth.

"""
is_smooth(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(i -> is_smooth(P, i), 1 : N)


@doc raw"""
    picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer}

The product of all local multiplicities of `P` divided by the global
multiplicity. For ldp polygons, this equals the index of the Picard group inside
the divisor class group of the associated toric surface, see [Spr24](@cite). 

"""
picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
prod([multiplicity(P,i) for i = 1 : N]) ÷ multiplicity(P)


@doc raw"""
    gorenstein_index(P :: RationalPolygon{T}, i :: Int)

The multiplicity of `P` divided by `gcd(w[2] - v[2], v[1] - w[1])`, where `v`
and `w` are the `i`-th and `i+1`-th scaled vertices of `P` respectively. For
ldp polygons, this equals the local gorenstein at the toric fixed point
associated to the `i`-th and `i+1`-th ray of `P`, see e.g. Lemma 3.9
of [HHHS22](@cite).

"""
function gorenstein_index(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
    v1, v2 = scaled_vertex(P,i), scaled_vertex(P,i+1)
    return det(v1,v2) ÷ gcd(v2[2] - v1[2], v1[1] - v2[1])
end


@doc raw"""
    gorenstein_index(P :: RationalPolygon{T}) where {T <: Integer}

The least common multiple of the local gorenstein indices of `P`. For ldp
polygons, this equals the gorenstein index of the associated toric surface.

"""
function gorenstein_index(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    g = 1
    for i = 1 : N
        g = lcm(g, gorenstein_index(P,i))
    end
    return g
end

function log_canonicities(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer}
    V = vertex_matrix(P)
    v1, v2 = V[:,mod(i,1:N)], V[:,mod(i+1,1:N)]
    return discrepancies(primitivize(v1), primitivize(v2))
end

log_canonicity(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} =
minimum(log_canonicities(P,i))


@doc raw"""
    log_canonicity(P :: RationalPolygon)

Given a `k`-rational polygon `P`, return the maximal rational number 0 < ϵ ≤ 1
such that ε*P contains only one `k`-rational point in its interior (the
origin). For an ldp polygon, this equals the maximal rational number 0 < ε ≤ 1
such that the associated toric surface is ε-log canonical.

"""
log_canonicity(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
minimum([log_canonicity(P,i) for i = 1 : N])


@doc raw"""
    toric_prime_divisor_self_intersection(P :: RationalPolygon, i :: Int)

Writing `u`, `v` and `w` for the `i-1`-th, `i`-th and `i+1`-th scaled vertex of
`P` respectively, return `det(w,u) // (det(u,v) * det(v,w))`. For ldp polygons,
this equals the self intersection number of the `i`-th toric prime divisor, see
e.g. Summary 3.2 of [HHS23](@cite).

"""
function toric_prime_divisor_self_intersection(P :: RationalPolygon{T,N}, i :: Int) where {N, T <: Integer}
    u, v, w = scaled_vertex(P,i-1), scaled_vertex(P,i), scaled_vertex(P,i+1)
    return det(w,u) // (det(u,v) * det(v,w))
end


@doc raw"""
    toric_prime_divisor_adjacent_intersection(P :: RationalPolygon, i :: Int)

Writing `v` and `w` for the `i`-th and `i+1`-th scaled vertex of `P`, return `1
// det(v,w)`. For ldp polygons, this equals the intersection number between the
`i`-th and `i+1`-th toric prime divisors.

"""
toric_prime_divisor_adjacent_intersection(P :: RationalPolygon{T,N}, i :: Int) where {N, T <: Integer} =
1 // det(scaled_vertex(P,i), scaled_vertex(P,i+1))


@doc raw"""
    degree(P :: RationalPolygon)

For ldp polygons, return the self intersection number of an anticanonocal
divisor of the associated toric surface.

# Example

The projective plane has degree 9.

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,0),(0,1),(-1,-1)])
Rational polygon of rationality 1 with 3 vertices.

julia> degree(P)
9//1
```

"""
degree(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
sum([toric_prime_divisor_self_intersection(P,i) + 2 * toric_prime_divisor_adjacent_intersection(P,i) for i = 1 : N])
