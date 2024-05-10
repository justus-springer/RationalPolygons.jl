
function _generic_k_rational_points(k :: T, P :: RationalPolygon{T,N}; mode :: Symbol = :interior) where {N,T <: Integer}

    mode ∈ [:interior, :integer_hull, :boundary, :all] || error("mode must be one of :interior, :integer_hull, :boundary and :all")

    ymin, ymax = minimum([v[2] for v ∈ P]), maximum([v[2] for v ∈ P])

    if mode == :interior
        lower_bound_y = floor_k_rational(k, ymin+1//k)
        upper_bound_y = ceil_k_rational(k, ymax-1//k)
    else
        lower_bound_y = ceil_k_rational(k, ymin)
        upper_bound_y = floor_k_rational(k, ymax)
    end

    res = RationalPoint{T}[]
    for y = lower_bound_y : 1 // k : upper_bound_y
        i0 = first(filter(i -> P[i][2] ≥ y ≥ P[i+1][2] && P[i][2] > P[i+1][2], 1 : N))
        incoming_line = line_through_points(P[i0], P[i0+1])
        incoming_point = intersection_point(horizontal_line(y), incoming_line)
        xmin = incoming_point[1]

        i1 = first(filter(i -> P[i][2] ≤ y ≤ P[i+1][2] && P[i][2] < P[i+1][2], 1 : N))
        outgoing_line = line_through_points(P[i1], P[i1+1])
        outgoing_point = intersection_point(horizontal_line(y), outgoing_line)
        xmax = outgoing_point[1]

        if mode == :interior
            lower_bound_x = floor_k_rational(k, xmin+1//k)
            upper_bound_x = ceil_k_rational(k, xmax-1//k)
        else
            lower_bound_x = ceil_k_rational(k, xmin)
            upper_bound_x = floor_k_rational(k, xmax)
        end

        if mode ∈ [:interior, :all]
            new_points = [RationalPoint{T}(x,y) for x = lower_bound_x : 1 // k : upper_bound_x]
        elseif mode ∈ [:boundary, :integer_hull]
            if lower_bound_x < upper_bound_x && (lower_bound_x, y) ∈ collect(P)
                new_points = [RationalPoint{T}(x,y) for x = lower_bound_x : 1 // k : upper_bound_x]
            elseif lower_bound_x < upper_bound_x
                new_points = [RationalPoint{T}(lower_bound_x, y), RationalPoint{T}(upper_bound_x, y)]
            elseif lower_bound_x == upper_bound_x
                new_points = [RationalPoint{T}(lower_bound_x, y)]
            else
                new_points = RationalPoint{T}[]
            end
        end

        append!(res, new_points)

    end

    return res

end


@doc raw"""
    boundary_k_rational_points(k :: T, P :: RationalPolygon{T,N}) where {T <: Integer} 

Return all `k`-rational points on the boundary of `P`.

"""
boundary_k_rational_points(k :: T, P :: RationalPolygon{T,N}) where {N,T <: Integer} =
unique(vcat([k_rational_points_on_line_segment(k, P[i], P[i+1]; interior = false) for i = 1 : N]...))


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
boundary_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(boundary_k_rational_points(T(1), P))


@doc raw"""
    number_of_boundary_lattice_points(P :: RationalPolygon)

Return the number of lattice points on the boundary of `P`.

"""
number_of_boundary_lattice_points(P :: RationalPolygon) =
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
interior_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(interior_k_rational_points(T(1), P))


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
lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
numerator.(k_rational_points(T(1), P))


@doc raw"""
    number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer}
    
Return the number of lattice points in `P`.

"""
number_of_lattice_points(P :: RationalPolygon{T}) where {T <: Integer} =
length(lattice_points(P))


k_rational_hull(k :: T, P :: RationalPolygon{T}) where {T <: Integer} =
convex_hull(_generic_k_rational_points(k, P; mode = :integer_hull), k)

integer_hull(P :: RationalPolygon{T}) where {T <: Integer} =
k_rational_hull(T(1), P)

