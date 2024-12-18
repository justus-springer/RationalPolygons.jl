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
    is_fano(P :: RationalPolygon)

Check whether `P` is a fano polygon, i.e. is primitive and contains the origin
in its interior.

"""
is_fano(P :: RationalPolygon) = contains_origin_in_interior(P) && is_primitive(P)


@doc raw"""
    dual(P :: RationalPolygon{T}) where {T <: Integer}

Return the dual of a polygon `P`. Note that `P` must contain the origin in its interior.

"""
function dual(P :: RationalPolygon{T}) where {T <: Integer}
    contains_origin_in_interior(P) || error("this polygon does not contain the origin in its interior")
    Hs = affine_halfplanes(P)
    return convex_hull([normal_vector(H) // translation(H) for H ∈ Hs])
end

local_class_group_order(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} =
det(scaled_vertex(P,i), scaled_vertex(P,i+1))


class_group_torsion_order(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
gcd([local_class_group_order(P,i) for i = 1 : N])


picard_index(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
prod([local_class_group_order(P,i) for i = 1 : N]) ÷ class_group_torsion_order(P)

@doc raw"""
    local_gorenstein_index(P :: RationalPolygon{T}, i :: Int)

Return the local gorenstein index at the `i`-th facet of `P`. This is onlywell defined for fano polygons.

"""
function local_gorenstein_index(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
    v1, v2 = scaled_vertex(P,i), scaled_vertex(P,i+1)
    return det(v1,v2) ÷ gcd(v2[2] - v1[2], v1[1] - v2[1])
end


@doc raw"""
    gorenstein_index(P :: RationalPolygon{T}) where {T <: Integer}

Return the gorenstein index of a `k`-rational polygon `P`, i.e. the smallest
positive integer `i` such that `i * dual(k * P)` is a lattice polygon.

"""
function gorenstein_index(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    g = 1
    for i = 1 : N
        g = lcm(g, local_gorenstein_index(P,i))
    end
    return g
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
origin). For a fano polygon, this equals the maximal rational number 0 < ε ≤ 1
such that the associated toric surface is ε-log canonical.

"""
log_canonicity(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
minimum([log_canonicity(P,i) for i = 1 : N])

function toric_prime_divisor_self_intersection(P :: RationalPolygon{T,N}, i :: Int) where {N, T <: Integer}
    u, v, w = scaled_vertex(P,i-1), scaled_vertex(P,i), scaled_vertex(P,i+1)
    return det(w,u) // (det(u,v) * det(v,w))
end

toric_prime_divisor_adjacent_intersection(P :: RationalPolygon{T,N}, i :: Int) where {N, T <: Integer} =
1 // det(scaled_vertex(P,i), scaled_vertex(P,i+1))

degree(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
sum([toric_prime_divisor_self_intersection(P,i) + 2 * toric_prime_divisor_adjacent_intersection(P,i) for i = 1 : N])
