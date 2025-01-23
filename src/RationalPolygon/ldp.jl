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
    picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer}

The product of all local multiplicities of `P` divided by the global
multiplicity. For ldp polygons, this equals the index of the Picard group inside
the divisor class group of the associated toric surface, see [Spr24](@cite). 

"""
picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
prod([multiplicity(P,i) for i = 1 : N]) ÷ multiplicity(P)


gorenstein_index(v :: LatticePoint{T}, w :: LatticePoint{T}) where {T <: Integer} = det(v,w) ÷ gcd(w[2] - v[2], v[1] - w[1])

@doc raw"""
    gorenstein_index(P :: RationalPolygon{T}, i :: Int)

The multiplicity of `P` divided by `gcd(w[2] - v[2], v[1] - w[1])`, where `v`
and `w` are the `i`-th and `i+1`-th scaled vertices of `P` respectively. For
ldp polygons, this equals the local gorenstein at the toric fixed point
associated to the `i`-th and `i+1`-th ray of `P`, see e.g. Lemma 3.9
of [HHHS22](@cite).

"""
gorenstein_index(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} =
gorenstein_index(scaled_vertex(P,i), scaled_vertex(P,i+1))


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

function degree_matrix(P :: RationalPolygon{T,N}) where {N, T <: Integer}
    S, U, _ = snf_with_transform(transpose(vertex_matrix(P)))
    d = S[2,2] # equals multiplicity(P)
    Q_free = SMatrix{N-2,N,T}(U[3:end,:])
    Q_torsion = SVector{N,T}(U[2,:] .% d)
    return Q_free, Q_torsion
end

degree_matrix_free_part(P :: RationalPolygon) = degree_matrix(P)[1]

degree_matrix_torsion_part(P :: RationalPolygon) = degree_matrix(P)[2]

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

function gorenstein_matrix(P :: RationalPolygon{T,N}) where {T <: Integer, N}
    A = gorenstein_coefficients(P)
    g = gorenstein_index(P)
    return SMatrix{N,N,T}([g - (mod(j-i+1,1:N) ≤ N-2 ? A[i,mod(j-i+1,1:N)] : 0) for i = 1:N, j = 1:N])
end

function log_canonicities(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer}
    V = vertex_matrix(P)
    v1, v2 = V[:,mod(i,1:N)], V[:,mod(i+1,1:N)]
    L = line_through_points(v1,v2)
    res = Rational{T}[]
    for p ∈ hilbert_basis(primitivize(v1),primitivize(v2))
        push!(res, norm_ratio(primitivize(p), intersection_point(L, line_through_points(zero(RationalPoint{T}), p))))
    end
    return res
end

log_canonicities(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
[log_canonicities(P,i) for i = 1 : N]

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
