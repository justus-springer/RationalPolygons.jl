
function _generic_k_rational_points(k :: T, P :: RationalPolygon{T}; mode :: Symbol = :interior) where {T <: Integer}
    vs, es = vertices(P), edges(P)

    mode ∈ [:interior, :integer_hull, :boundary, :all] || error("mode must be one of :interior, :integer_hull, :boundary and :all")

    ymin, ymax = minimum([v[2] for v ∈ vs]), maximum([v[2] for v ∈ vs])

    if mode == :interior
        lower_bound_y = floor_k_rational(k, ymin+1//k)
        upper_bound_y = ceil_k_rational(k, ymax-1//k)
    else
        lower_bound_y = ceil_k_rational(k, ymin)
        upper_bound_y = floor_k_rational(k, ymax)
    end

    res = RationalPoint{T}[]
    for y = lower_bound_y : 1 // k : upper_bound_y
        incoming_edge = first(filter(e -> e[1][2] ≥ y ≥ e[2][2] && e[1][2] > e[2][2], es))
        incoming_line = LineThroughPoints(incoming_edge[1], incoming_edge[2])
        incoming_point = intersection_point(HorizontalLine(y), incoming_line)
        xmin = incoming_point[1]

        outgoing_edge = first(filter(e -> e[1][2] ≤ y ≤ e[2][2] && e[1][2] < e[2][2], es))
        outgoing_line = LineThroughPoints(outgoing_edge[1], outgoing_edge[2])
        outgoing_point = intersection_point(HorizontalLine(y), outgoing_line)
        xmax = outgoing_point[1]

        if mode == :interior
            lower_bound_x = floor_k_rational(k, xmin+1//k)
            upper_bound_x = ceil_k_rational(k, xmax-1//k)
        else
            lower_bound_x = ceil_k_rational(k, xmin)
            upper_bound_x = floor_k_rational(k, xmax)
        end

        if mode ∈ [:interior, :all]
            new_points = [(x,y) for x = lower_bound_x : 1 // k : upper_bound_x]
        elseif mode ∈ [:boundary, :integer_hull]
            if lower_bound_x < upper_bound_x && (lower_bound_x, y) ∈ vertices(P)
                new_points = [(x,y) for x = lower_bound_x : 1 // k : upper_bound_x]
            elseif lower_bound_x < upper_bound_x
                new_points = [(lower_bound_x, y), (upper_bound_x, y)]
            elseif lower_bound_x == upper_bound_x
                new_points = [(lower_bound_x, y)]
            else
                new_points = RationalPoint{T}[]
            end
        end

        append!(res, new_points)

    end

    return res

end


@doc raw"""
    boundary_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} 

Return all `k`-rational points on the boundary of `P`.

"""
boundary_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
unique(vcat([k_rational_points_on_line_segment(k, e[1], e[2]; interior = false) for e ∈ edges(P)]...))


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
@attr boundary_lattice_points(P :: RationalPolygon) =
numerator.(boundary_k_rational_points(1, P))


@doc raw"""
    number_of_boundary_lattice_points(P :: RationalPolygon)

Return the number of lattice points on the boundary of `P`.

"""
@attr number_of_boundary_lattice_points(P :: RationalPolygon) =
length(boundary_lattice_points(P))


@doc raw"""
    interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return all `k`-rational points in the interior of `P`.

"""
interior_k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
_generic_k_rational_points(k, P; mode = :interior)


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
@attr interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(interior_k_rational_points(1, P))


@doc raw"""
    number_of_interior_lattice_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}

Return the number of `k`-rational points in the interior of `P`.

"""
@attr number_of_interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
length(interior_lattice_points(P))


@doc raw"""
    k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer}
    
Return all `k`-rational points in `P`.

"""
k_rational_points(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
_generic_k_rational_points(k, P; mode = :all)


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
@attr lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(k_rational_points(1, P))


@doc raw"""
    number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}
    
Return the number of lattice points in `P`.

"""
@attr number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
length(lattice_points(P))


k_rational_hull(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
convex_hull(_generic_k_rational_points(k, P; mode = :integer_hull); rationality = k)

@attr integer_hull(P :: RationalPolygon{T}) where {T <: Integer} =
k_rational_hull(1, P)

