
@doc raw"""
    affine_halfplane(P :: RationalPolygon, i :: Int)

Return the `i`-th describing halfplane of `P`.

"""
affine_halfplane(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} =
affine_halfplane(line_through_points(P[i], P[i+1]))


@doc raw"""
    affine_halfplanes(P :: RationalPolygon)

Return the describing halfplanes of `P`.

"""
affine_halfplanes(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
[affine_halfplane(P, i) for i = 1 : N]


@doc raw"""
    Base.in(x :: Point{T}, P :: RationalPolygon{T}) where {T <: Integer}

Check whether a point `x` is contained in the polygon `P`.

"""
Base.in(x :: Point{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(H -> x ∈ H, affine_halfplanes(P))


@doc raw"""
    contains_in_interior(x :: Point{T}, P :: RationalPolygon{T}) where {T <: Integer}

Check whether a point `x` is contained in the interior of `P`.

"""
contains_in_interior(x :: Point{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(H -> contains_in_interior(x, H), affine_halfplanes(P))


@doc raw"""
    dim(P :: RationalPolygon)

Return the dimension of `P`. For empty polygons, this returns -1. Otherwise,
it returns 0, 1 or 2.

"""
dim(P :: RationalPolygon{T,0}) where {T <: Integer} = -1
dim(P :: RationalPolygon{T,1}) where {T <: Integer} = 0 
dim(P :: RationalPolygon{T,2}) where {T <: Integer} = 1 
dim(P :: RationalPolygon{T,N}) where {N,T <: Integer}  = 2 


@doc raw"""
    dimension_of_interior_integer_hull(P :: RationalPolygon)   

Return the dimension of the convex hull of the interior lattice points of `P`.

"""
dimension_of_interior_integer_hull(P :: RationalPolygon) =
dim(interior_integer_hull(P))


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
    normalized_area(P :: RationalPolygon)

Return the normalized area of a `k`-rational polygon. The result is an
integer, counting the number of standard `k`-rational triangles contained in
`P`.

"""
function normalized_area(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    V = vertex_matrix(P)
    return sum([abs(det(V[:,i+1] - V[:,1], V[:,i] - V[:,1])) for i = 2 : N -1])
end


@doc raw"""
    euclidian_area(P :: RationalPolygon)

Return the euclidian area of a rational polygon.

"""
euclidian_area(P :: RationalPolygon) = normalized_area(P) // (2 * rationality(P)^2)


@doc raw"""
    is_maximal(P :: RationalPolygon)

Check whether a rational polygon is maximal among all polygons sharing
the same rationality and number of interior lattice points.

"""
function is_maximal(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    k = rationality(P)
    Hs = affine_halfplanes(P)
    Q = intersect_halfplanes(Hs .- 1 // k)

    nonempty_edges_indices = filter(i -> !isempty(integral_points_on_line_segment(P[i], P[i+1])), 1 : N)
    integral_vertices_indices = filter(i -> is_integral(P[i]), 1 : N)

    for p ∈ boundary_k_rational_points(Q, k)
        if all(i -> p ∈ Hs[i], nonempty_edges_indices) &&
           all(i -> p ∈ Hs[mod(i-1,1:N)] || p ∈ Hs[i], integral_vertices_indices)
            return false
        end
    end

    return true

end


@doc raw"""
    move_out_edges(P :: RationalPolygon)

Given a `k`-rational polygon `P`, return the polygon obtained by moving out all
edges by `1 // k`.

"""
function move_out_edges(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    k = rationality(P)
    Hs = affine_halfplanes(P)
    return intersect_halfplanes(Hs .- one(T) // k)
end


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
    gorenstein_index(P :: RationalPolygon{T}) where {T <: Integer}

Return the gorenstein index of a `k`-rational polygon `P`, i.e. the smallest
positive integer `i` such that `i * dual(k * P)` is a lattice polygon.

"""
function gorenstein_index(P :: RationalPolygon{T}) where {T <: Integer}
    is_fano(P) || error("the gorenstein index is only defined for fano polygons")
    k = rationality(P)
    local_gorenstein_indices = [numerator(k*translation(H)) for H ∈ affine_halfplanes(P)]
    return lcm(local_gorenstein_indices)
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
