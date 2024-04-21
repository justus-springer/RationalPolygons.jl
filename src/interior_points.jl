@doc raw"""
    boundary_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} 

Return all `k`-rational points on the boundary of `P`.

"""
boundary_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
unique(vcat([k_rational_points_on_line_segment(k, e[1], e[2]; inclusive = true) for e ∈ edges(P)]...))


@doc raw"""
    number_of_boundary_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}    

Return the number of `k`-rational points on the boundary of `P`.

"""
number_of_boundary_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
length(boundary_k_rational_points(k,P))


@doc raw"""
    boundary_lattice_points(P :: RationalPolygon)

Return all lattice points on the boundary of `P`.

"""
boundary_lattice_points(P :: RationalPolygon) =
numerator.(boundary_k_rational_points(1, P))


@doc raw"""
    number_of_boundary_lattice_points(P :: RationalPolygon)

Return the number of lattice points on the boundary of `P`.

"""
number_of_boundary_lattice_points(P :: RationalPolygon) =
length(boundary_lattice_points(P))


function _generic_k_rational_points(k :: T, P :: RationalPolygon{T}; interior = true) where {T <: Integer}
    vs, es = vertices(P), edges(P)

    ymin, ymax = minimum([v[2] for v ∈ vs]), maximum([v[2] for v ∈ vs])
    lower_bound = interior ? floor_k_rational(k, ymin+1//k) : ceil_k_rational(k, ymin)
    upper_bound = interior ? ceil_k_rational(k, ymax-1//k) : floor_k_rational(k, ymax)

    res = RationalPoint[]
    for t = lower_bound : 1 // k : upper_bound
        incoming_edge = first(filter(e -> e[1][2] ≥ t ≥ e[2][2] && e[1][2] > e[2][2], es))
        incoming_line = LineThroughPoints(incoming_edge[1], incoming_edge[2])
        incoming_point = intersection_point(HorizontalLine(t), incoming_line)

        outgoing_edge = first(filter(e -> e[1][2] ≤ t ≤ e[2][2] && e[1][2] < e[2][2], es))
        outgoing_line = LineThroughPoints(outgoing_edge[1], outgoing_edge[2])
        outgoing_point = intersection_point(HorizontalLine(t), outgoing_line)

        append!(res, k_rational_points_on_line_segment(k, incoming_point, outgoing_point; interior))

    end

    return res

end


@doc raw"""
    interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return all `k`-rational points in the interior of `P`.

"""
interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
_generic_k_rational_points(k, P; interior = true)


@doc raw"""
    number_of_interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return the number of `k`-rational points in the interior of `P`.

"""
number_of_interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
length(interior_k_rational_points(k,P))


@doc raw"""
    interior_lattice_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return all `k`-rational points in the interior of `P`.

"""
interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(interior_k_rational_points(1, P))


@doc raw"""
    number_of_interior_lattice_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return the number of `k`-rational points in the interior of `P`.

"""
number_of_interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
length(interior_lattice_points(P))


@doc raw"""
    k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}
    
Return all `k`-rational points in `P`.

"""
k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
_generic_k_rational_points(k, P; interior = false)


@doc raw"""
    number_of_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return the number of `k`-rational points in `P`.

"""
number_of_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
length(k_rational_points(k,P))


@doc raw"""
    lattice_points(P :: RationalPolygon{T}) where {T <: Integer}
    
Return all lattice points in `P`.

"""
lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
k_rational_points(1, P)


@doc raw"""
    number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}
    
Return the number of lattice points in `P`.

"""
number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
length(lattice_points(P))
