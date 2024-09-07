
@doc raw"""
    generic_lattice_points(P :: RationalPolygon{T,N}, k :: T)

Return the lattice points in the `k`-fold of the rational polygon `P`. If the
keyword argument `interior = true` is passed, only the interior lattice points
are computed. If `only_count = true` is passed, only the total number of
lattice points is returned.
        
"""
function generic_lattice_points(
        P :: RationalPolygon{T,N},
        k :: T;
        interior :: Bool = false,
        only_count :: Bool = false) where {N,T <: Integer}

    vs = k .* vertices(P)
    ymin, ymax = minimum([v[2] for v ∈ vs]), maximum([v[2] for v ∈ vs])

    lby = interior ? floor(T,ymin+1) : ceil(T,ymin)
    uby = interior ? ceil(T,ymax-1) : floor(T,ymax)

    count = 0
    points = LatticePoint{T}[]
    for y = lby : uby
        i0 = first(filter(i -> vs[i][2] ≥ y ≥ vs[mod(i+1,1:N)][2] && 
                               vs[i][2] > vs[mod(i+1,1:N)][2], 1 : N))
        incoming_line = line_through_points(vs[i0], vs[mod(i0+1,1:N)])
        incoming_point = intersection_point(horizontal_line(y), incoming_line)
        xmin = incoming_point[1]

        i1 = first(filter(i -> vs[i][2] ≤ y ≤ vs[mod(i+1,1:N)][2] && vs[i][2] < vs[mod(i+1,1:N)][2], 1 : N))
        outgoing_line = line_through_points(vs[i1], vs[mod(i1+1,1:N)])
        outgoing_point = intersection_point(horizontal_line(y), outgoing_line)
        xmax = outgoing_point[1]

        if xmax < xmin
            xmax,xmin = xmin,xmax
        end

        lbx = interior ? floor(T,xmin+1) : ceil(T,xmin)
        ubx = interior ? ceil(T,xmax-1) : floor(T,xmax)
        count += ubx - lbx + 1
        if !only_count
            append!(points, [LatticePoint{T}(x,y) for x = lbx : ubx])
        end

    end

    return only_count ? count : points

end


@doc raw"""
    boundary_k_rational_points(P :: RationalPolygon{T,N}, k :: T) where {N,T <: Integer}

Return all `k`-rational points on the boundary of `P`.

"""
function boundary_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer}
    res = RationalPoint{T}[]
    for i = 1 : number_of_vertices(P)
        append!(res, k_rational_points_on_line_segment(k, P[i], P[i+1]; interior = false))
    end
    return unique(res)
end


@doc raw"""
    number_of_boundary_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer}

Return the number of `k`-rational points on the boundary of `P`.

"""
number_of_boundary_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer} =
length(boundary_k_rational_points(P,k))


@doc raw"""
    boundary_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}

Return all lattice points on the boundary of `P`.

"""
boundary_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(boundary_k_rational_points(P, one(T)))


@doc raw"""
    number_of_boundary_lattice_points(P :: RationalPolygon)

Return the number of lattice points on the boundary of `P`.

"""
number_of_boundary_lattice_points(P :: RationalPolygon) =
length(boundary_lattice_points(P))


@doc raw"""
    interior_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer}

Return all `k`-rational points in the interior of `P`.

"""
interior_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer} =
generic_lattice_points(P, k; interior = true) .// k


@doc raw"""
    number_of_interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return the number of `k`-rational points in the interior of `P`.

"""
number_of_interior_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer} =
generic_lattice_points(P, k; interior = true, only_count = true)


@doc raw"""
    interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}

Return all lattice points in the interior of `P`.

"""
interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
generic_lattice_points(P, one(T); interior = true)


@doc raw"""
    number_of_interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}

Return the number of lattice points in the interior of `P`.

"""
number_of_interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
generic_lattice_points(P, one(T), interior = true, only_count = true)


@doc raw"""
    k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}
    
Return all `k`-rational points in `P`.

"""
k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer} =
generic_lattice_points(P, k) .// k


@doc raw"""
    number_of_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer}

Return the number of `k`-rational points in `P`.

"""
number_of_k_rational_points(P :: RationalPolygon{T}, k :: T) where {T <: Integer} =
generic_lattice_points(P, k, only_count = true)


@doc raw"""
    lattice_points(P :: RationalPolygon{T}) where {T <: Integer}
    
Return all lattice points in `P`.

"""
lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
generic_lattice_points(P, one(T))


@doc raw"""
    number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}
    
Return the number of lattice points in `P`.

"""
number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
generic_lattice_points(P, one(T), only_count = true)


@doc raw"""
    k_rational_hull(P :: RationalPolygon{T}, k :: T) where {T <: Integer}

Return the convex hull of all `k`-rational points contained in `P`. If
`primitive = true` is passed, only the primitive `k`-rational points
are taken.

"""
function k_rational_hull(P :: RationalPolygon{T}, k :: T; primitive = false) where {T <: Integer}
    ps = k_rational_points(P, k)
    primitive && filter!(p -> is_primitive(k*p), ps)
    return convex_hull(ps, k)
end


@doc raw"""
    interior_k_rational_hull(P :: RationalPolygon{T}, k :: T) where {T <: Integer}

Return the convex hull of all interior `k`-rational points of `P`. If
`primitive = true` is passed, only the primitive `k`-rational points
are taken.

"""
function interior_k_rational_hull(P :: RationalPolygon{T}, k :: T; primitive = false) where {T <: Integer}
    ps = interior_k_rational_points(P, k)
    primitive && filter!(p -> is_primitive(k*p), ps)
    return convex_hull(ps, k)
end


@doc raw"""
    integer_hull(P :: RationalPolygon{T}) where {T <: Integer}

Return the convex hull of all lattice points in the interior of `P`. If
`primitive = true` is passed, only the primitive interior lattice points are
taken.

"""
integer_hull(P :: RationalPolygon{T}; primitive = false) where {T <: Integer} =
k_rational_hull(P, one(T); primitive)


@doc raw"""
    interior_integer_hull(P :: RationalPolygon{T}) where {T <: Integer}

Return the convex hull of all interior lattice points of `P`. If `primitive =
true` is passed, only the primitive lattice points are taken.

"""
interior_integer_hull(P :: RationalPolygon{T}; primitive = false) where {T <: Integer} =
interior_k_rational_hull(P, one(T); primitive)
