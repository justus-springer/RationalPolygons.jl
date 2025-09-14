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

Check whether `P` is LDP, i.e. is primitive and contains the origin
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

Return the dual of a polygon `P`. Throws an error if `P` does not contain the origin in its interior.

"""
function dual(P :: RationalPolygon{T}) where {T <: Integer}
    contains_origin_in_interior(P) || error("this polygon does not contain the origin in its interior")
    Hs = affine_halfplanes(P)
    return convex_hull([normal_vector(H) // translation(H) for H ∈ Hs])
end


@doc raw"""
    multiplicity(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

The index of the sublattice spanned by the `i`-th and `i+1`-th scaled vertex of
`P` (i.e. the determinant of those two vertices). For LDP polygons, this equals
the order of the local class group associated with the toric fixed point
associated to the `i`-th and `i+1`-th ray.

"""
multiplicity(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} =
det(scaled_vertex(P,i), scaled_vertex(P,i+1))


@doc raw"""
    multiplicity(P :: RationalPolygon{T,N}) where {N,T <: Integer}

The order of the sublattice spanned by the scaled vertices of `P`. For LDP
polygons, this equals the order of the torsion part of the divisor class group
of the associated toric surface.

"""
multiplicity(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
gcd([multiplicity(P,i) for i = 1 : N])


@doc raw"""
    is_smooth(P :: RationalPolygon, i :: Int)

Check whether the cone spanned by the `i`-th and `i+1`-th vertex of `P` is
regular, i.e. generates the entire lattice. For LDP polygons, this means that
the toric fixed point associated to the `i`-th and `i+1`-th ray is smooth.

"""
is_smooth(P :: RationalPolygon, i :: Int) =
multiplicity(P, i) == 1


@doc raw"""
    is_smooth(P :: RationalPolygon)

Check whether all cones of the face fan of `P` are regular. For LDP polygons,
this means that the associated toric surface is smooth.

"""
is_smooth(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(i -> is_smooth(P, i), 1 : N)


@doc raw"""
    picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer}

The product of all local multiplicities of `P` divided by the global
multiplicity. For LDP polygons, this equals the index of the Picard group inside
the divisor class group of the associated toric surface, see Proposition ``\ref{prp:picard_index_formula_ldp_polygons}``

"""
picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
prod([multiplicity(P,i) for i = 1 : N]) ÷ multiplicity(P)


gorenstein_index(v :: LatticePoint{T}, w :: LatticePoint{T}) where {T <: Integer} = det(v,w) ÷ gcd(w[2] - v[2], v[1] - w[1])

@doc raw"""
    gorenstein_index(P :: RationalPolygon{T}, i :: Int)

The multiplicity of `P` divided by `gcd(w[2] - v[2], v[1] - w[1])`, where `v`
and `w` are the `i`-th and `i+1`-th scaled vertices of `P` respectively. For
LDP polygons, this equals the local Gorenstein at the toric fixed point
associated to the `i`-th and `i+1`-th ray of `P`, see Proposition ``\ref{prp:gorenstein_index_formula}``

"""
gorenstein_index(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} =
gorenstein_index(scaled_vertex(P,i), scaled_vertex(P,i+1))


@doc raw"""
    gorenstein_index(P :: RationalPolygon{T}) where {T <: Integer}

The least common multiple of the local Gorenstein indices of `P`. For LDP
polygons, this equals the Gorenstein index of the associated toric surface.

"""
function gorenstein_index(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    g = 1
    for i = 1 : N
        g = lcm(g, gorenstein_index(P,i))
    end
    return g
end

@doc raw"""
    degree_matrix(P :: RationalPolygon{T,N}) where {N, T <: Integer}

Return a tuple ``(Q_0, Q_1)`` where ``Q_0`` is the free part and ``Q_1``
the torsion part of the degree matrix associated to ``P``.

"""
function degree_matrix(P :: RationalPolygon{T,N}) where {N, T <: Integer}
    S, U, _ = snf_with_transform(transpose(vertex_matrix(P)))
    d = S[2,2] # equals multiplicity(P)
    Q_free = SMatrix{N-2,N,T}(U[3:end,:])
    Q_torsion = SVector{N,T}([mod(U[2,i], 0:d-1) for i = 1 : N])
    return Q_free, Q_torsion
end


@doc raw"""
    degree_matrix_free_part(P :: RationalPolygon)

Return the free part of the degree matrix associated to ``P``.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,0), (2,5), (-4,-5), (-1,-5)])
Rational polygon of rationality 1 with 4 vertices.

julia> degree_matrix_free_part(P)
2×4 StaticArraysCore.SMatrix{2, 4, Int64, 8} with indices SOneTo(2)×SOneTo(4):
 1  -1  3  0
 1   0  2  1
```

"""
degree_matrix_free_part(P :: RationalPolygon) = degree_matrix(P)[1]


@doc raw"""
    degree_matrix_torsion_part(P :: RationalPolygon)

Return the torsion part of the degree matrix associated to ``P``.
If ``P`` has ``N`` vertices and has multiplicity ``\mu``, the result is a static vector of length ``N`` whose 
entries are all between ``0`` and ``\mu-1``.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,0), (2,5), (-4,-5), (-1,-5)])
Rational polygon of rationality 1 with 4 vertices.

julia> multiplicity(P)
5

julia> degree_matrix_torsion_part(P)
4-element StaticArraysCore.SVector{4, Int64} with indices SOneTo(4):
 4
 0
 1
 0
```

"""
degree_matrix_torsion_part(P :: RationalPolygon) = degree_matrix(P)[2]


@doc raw"""
    gorenstein_coefficients(P :: RationalPolygon{T,N})

Return the Gorenstein coefficients of an LDP polygon with ``N`` vertices.
This is an integral matrix ``A = (a_{ij}) \in \ZZ^{N \times (N-2)}`` such that
``\iota w = \sum_{j=1}^{N-2} a_{ij} w_{j+i-1}``, where ``\iota`` is the
Gorenstein index, ``w_i`` are the columns of the free part of the degree matrix,
and ``w = w_1 + \dots + w_N`` is the class of the anticanonical divisor. See
Definition ``\ref{def:gorenstein_coefficients}``.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,0), (2,5), (-4,-5), (-1,-5)])
Rational polygon of rationality 1 with 4 vertices.

julia> gorenstein_coefficients(P)
4×2 StaticArraysCore.SMatrix{4, 2, Int64, 8} with indices SOneTo(4)×SOneTo(2):
 20   5
 15  10
  5  10
  5  15
```


"""
function gorenstein_coefficients(P :: RationalPolygon{T,N}) where {T <: Integer, N}
    Q = degree_matrix_free_part(P)
    g = gorenstein_index(P)
    w = g * sum([Q[:, i] for i = 1 : N])
    As = SVector{N-2,T}[]
    for i = 1 : N
        js = map(j -> mod(j, 1:N), i:i+N-3)
        M = SMatrix{N-2,N-2,T}(Q[:, js])
        push!(As, solve_unique_integer_solution(M, w))
    end
    return transpose(hcat(As...))
end

gorenstein_matrix(ι :: T, A :: SMatrix{N,M,T}) where {T <: Integer, N, M} =
SMatrix{N,N,T}([(mod(j-i+1,1:N) ≤ N-2 ? ι - A[i,mod(j-i+1,1:N)] : ι) for i = 1:N, j = 1:N])

@doc raw"""
    gorenstein_matrix(P :: RationalPolygon{T,N}) where {T <: Integer, N}

Return the Gorenstein matrix of an LDP polygon. See Definition ``\ref{def:gorenstein_matrix}``.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,0), (2,5), (-4,-5), (-1,-5)])
Rational polygon of rationality 1 with 4 vertices.

julia> gorenstein_matrix(P)
4×4 StaticArraysCore.SMatrix{4, 4, Int64, 16} with indices SOneTo(4)×SOneTo(4):
 -15    0   5   5
   5  -10  -5   5
   5    5   0  -5
 -10    5   5   0
```

"""
gorenstein_matrix(P :: RationalPolygon{T,N}) where {T <: Integer, N} =
gorenstein_matrix(gorenstein_index(P), gorenstein_coefficients(P))

function log_canonicities(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer}
    V = vertex_matrix(P)
    v1, v2 = V[:,mod(i,1:N)], V[:,mod(i+1,1:N)]
    return discrepancies(primitivize(v1), primitivize(v2))
end

log_canonicity(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} =
minimum(log_canonicities(P,i))


@doc raw"""
    log_canonicity(P :: RationalPolygon)

Given a ``k``-rational polygon ``P``, return the maximal rational number ``0 < \varepsilon \leq 1``
such that ``\varepsilon*P`` contains only one ``k``-rational point in its interior (the
origin). For an LDP polygon, this equals the maximal rational number ``0 < \varepsilon \leq 1``
such that the associated toric surface is ``\varepsilon``-log canonical.

"""
log_canonicity(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
minimum([log_canonicity(P,i) for i = 1 : N])


@doc raw"""
    toric_prime_divisor_self_intersection(P :: RationalPolygon, i :: Int)

Writing `u`, `v` and `w` for the `i-1`-th, `i`-th and `i+1`-th scaled vertex of
`P` respectively, return `det(w,u) // (det(u,v) * det(v,w))`. For LDP polygons,
this equals the self intersection number of the `i`-th toric prime divisor, see
e.g. Summary 3.2 of [HaHaSp25](@cite).

"""
function toric_prime_divisor_self_intersection(P :: RationalPolygon{T,N}, i :: Int) where {N, T <: Integer}
    u, v, w = scaled_vertex(P,i-1), scaled_vertex(P,i), scaled_vertex(P,i+1)
    return det(w,u) // (det(u,v) * det(v,w))
end


@doc raw"""
    toric_prime_divisor_adjacent_intersection(P :: RationalPolygon, i :: Int)

Writing `v` and `w` for the `i`-th and `i+1`-th scaled vertex of `P`, return `1
// det(v,w)`. For LDP polygons, this equals the intersection number between the
`i`-th and `i+1`-th toric prime divisors.

"""
toric_prime_divisor_adjacent_intersection(P :: RationalPolygon{T,N}, i :: Int) where {N, T <: Integer} =
1 // det(scaled_vertex(P,i), scaled_vertex(P,i+1))


@doc raw"""
    degree(P :: RationalPolygon)

For LDP polygons, return the self intersection number of an anticanonocal
divisor of the associated toric surface.

# Example:

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
