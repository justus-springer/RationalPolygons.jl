
@doc raw"""
    width(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}

Return the lattice width of `P` in direction `w`.

"""
function width(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}
    values_on_vertices = [dot(v,w) for v ∈ vertices(P)]
    return maximum(values_on_vertices) - minimum(values_on_vertices)
end


@doc raw"""
    width_with_direction_vectors(P :: RationalPolygon{T}) where {T <: Integer}

Return the lattice width of `P` together with the list of direction vectors
that realize this width.

"""
function width_with_direction_vectors(P :: RationalPolygon{T}) where {T <: Integer}
    c = min(width(P, RationalPoint{T}(0,1)), width(P, RationalPoint{T}(1,0)))
    vs = lattice_points(c * dual(P - P))
    filter!(v -> v[1] > 0 || (v[1] == 0 && v[2] > 0), vs)
    widths = [width(P, v) for v ∈ vs]
    w = minimum(widths)
    direction_vectors = [vs[i] for i = 1 : length(vs) if widths[i] == w]
    return (w,direction_vectors)
end


@doc raw"""
    width(P :: RationalPolygon)

Return the lattice width of `P`.

"""
width(P :: RationalPolygon) = width_with_direction_vectors(P)[1]


@doc raw"""
    scaled_width(P :: RationalPolygon)

Return the width of the scaled lattice polygon `k * P`, where `k` is the
rationality of `P`. This equals `k * width(P)` and is always an integer.

"""
scaled_width(P :: RationalPolygon) = numerator(rationality(P) * width(P))


@doc raw"""
    width_direction_vectors(P :: RationalPolygon)

Return the lattice width direction vectors of `P`, i.e. those directions that
realize the lattice width of `P`.

"""
width_direction_vectors(P :: RationalPolygon) = width_with_direction_vectors(P)[2]

function adjust_to_width_direction(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}
    is_primitive(w) || error("the given direction vector is not primitive")
    k = rationality(P)
    _, x, y = gcdx(w[1],w[2])
    U = Matrix2{T}(w[1],-y,w[2],x)
    Q = U * P
    xmin = minimum(vertex_matrix(Q)[1,:])
    b = LatticePoint(-k * fld(xmin,k), 0)
    return Q + b
end

struct LatticeWidthData{T <: Integer}
    interval_of_nonzero_vertical_slice_length :: Tuple{Rational{T},Rational{T}}
    interval_of_longest_vertical_slice_length :: Tuple{Rational{T},Rational{T}}
end

number_of_interior_integral_vertical_lines(lwd :: LatticeWidthData{T}) where {T <: Integer} =
ceil(T, lwd.interval_of_nonzero_vertical_slice_length[2]) - 
floor(T, lwd.interval_of_nonzero_vertical_slice_length[1]) - 1

position_of_longest_vertical_slice_length(lwd :: LatticeWidthData) =
(lwd.interval_of_longest_vertical_slice_length[1] +
 lwd.interval_of_longest_vertical_slice_length[2]) // 2

function lattice_width_data(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}
    Q = adjust_to_width_direction(P, w)
    xs = [v[1] for v ∈ vertices(Q)]
    sort!(unique!(xs))
    vertical_slice_lengths = [slice_length(vertical_line(x), Q) for x ∈ xs]
    longest_vertical_slice_length = maximum(vertical_slice_lengths)
    max_positions = [xs[i] for i = 1 : length(xs) if vertical_slice_lengths[i] == longest_vertical_slice_length]

    return LatticeWidthData{T}((minimum(xs), maximum(xs)), 
                               (minimum(max_positions), maximum(max_positions)))

end

number_of_interior_integral_vertical_lines(P :: RationalPolygon, w :: Point{T}) where {T <: Integer} =
number_of_interior_integral_vertical_lines(lattice_width_data(P, w))

position_of_longest_vertical_slice_length(P :: RationalPolygon, w :: Point{T}) where {T <: Integer} =
position_of_longest_vertical_slice_length(lattice_width_data(P, w))

lattice_width_datas(P :: RationalPolygon{T}) where {T <: Integer} =
lattice_width_data.(P, width_direction_vectors(P))

numbers_of_interior_integral_vertical_lines(P :: RationalPolygon{T}) where {T <: Integer} =
[number_of_interior_integral_vertical_lines(P, w) for w ∈ width_direction_vectors(P)]

positions_of_longest_vertical_slice_length(P :: RationalPolygon{T}) where {T <: Integer} =
[position_of_longest_vertical_slice_length(P, w) for w ∈ width_direction_vectors(P)]
