
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
    normalized_area(P :: RationalPolygon)

Return the normalized area of a `k`-rational polygon, i.e. `2k^2` times its
euclidian area. The result is always an integer, counting the number of
standard `k`-rational triangles contained in `P`.

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
