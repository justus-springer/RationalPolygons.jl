
@doc raw"""
    width(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}

Return the lattice width of `P` in direction `w`.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,1),(1,-2),(-4,2),(-2,2)],2)
Rational polygon of rationality 2 with 4 vertices.

julia> width(P, Point(0,1))
2//1

julia> width(P, Point(1,0))
5//2
```

"""
function width(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}
    values_on_vertices = [dot(v,w) for v ∈ vertices(P)]
    return maximum(values_on_vertices) - minimum(values_on_vertices)
end


@doc raw"""
    all_direction_vectors_with_width_less_than(P :: RationalPolygon{T}, c :: Rational{T}) where {T <: Integer}

Return all direction vectors in which the width of `P` is less than or equal to a given constant.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,1),(1,-2),(-4,2),(-2,2)],2)
Rational polygon of rationality 2 with 4 vertices.

julia> all_direction_vectors_with_width_less_than(P, 3//1)
4-element Vector{StaticArraysCore.SVector{2, Int64}}:
 [1, 0]
 [0, 1]
 [1, 1]
 [1, 2]
```

"""
function all_direction_vectors_with_width_less_than(P :: RationalPolygon{T}, c :: Rational{T}) where {T <: Integer}
    vs = lattice_points(c * dual(P - P))
    filter!(v -> v[1] > 0 || (v[1] == 0 && v[2] > 0), vs)
    return vs
end


@doc raw"""
    width_with_direction_vectors(P :: RationalPolygon{T}) where {T <: Integer}

Return the lattice width of `P` together with the list of direction vectors
that realize this width.


"""
function width_with_direction_vectors(P :: RationalPolygon{T}) where {T <: Integer}
    c = min(width(P, RationalPoint{T}(0,1)), width(P, RationalPoint{T}(1,0)))
    vs = all_direction_vectors_with_width_less_than(P, c)
    widths = [width(P, v) for v ∈ vs]
    w = minimum(widths)
    direction_vectors = [vs[i] for i = 1 : length(vs) if widths[i] == w]
    return (w,direction_vectors)
end


@doc raw"""
    width(P :: RationalPolygon)

Return the lattice width of `P`.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,1),(1,-2),(-4,2),(-2,2)],2)
Rational polygon of rationality 2 with 4 vertices.

julia> width(P)
2//1
```

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

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,1),(1,-2),(-4,2),(-2,2)],2)
Rational polygon of rationality 2 with 4 vertices.

julia> width_direction_vectors(P)
2-element Vector{StaticArraysCore.SVector{2, Int64}}:
 [0, 1]
 [1, 1]
```

"""
width_direction_vectors(P :: RationalPolygon) = width_with_direction_vectors(P)[2]


@doc raw"""
    adjust_to_width_direction(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}

Apply an affine unimodular transformation to `P` that transforms the given
width direction vector to (1,0), see Lemma 2.10 of [Boh23](@cite).

"""
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


@doc raw"""
    number_of_interior_integral_lines(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}

Return the number of interior integral lines of `P` with respect to a given
direction vector `w`.

# Example

```jldoctest
julia> P = convex_hull(RationalPoint{Int}[(1//3,-1),(4//3,2),(2//3,2),(-4//3,-1)])
Rational polygon of rationality 3 with 4 vertices.

julia> number_of_interior_integral_lines(P, LatticePoint(1,0))
3

julia> number_of_interior_integral_lines(P, LatticePoint(0,1))
2
```

"""
function number_of_interior_integral_lines(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}
    k = rationality(P)
    Q = adjust_to_width_direction(P, w)
    xmin = minimum(vertex_matrix(Q)[1,:]) // k
    xmax = maximum(vertex_matrix(Q)[1,:]) // k
    return ceil(T,xmax) - floor(T,xmin) - 1
end


@doc raw"""
    minimal_number_of_interior_integral_lines(P :: RationalPolygon{T}) where {T <: Integer}

Return the minimal number of interior integral lines of `P`.

# Example

```jldoctest
julia> P = convex_hull(RationalPoint{Int}[(1//3,-1),(4//3,2),(2//3,2),(-4//3,-1)])
Rational polygon of rationality 3 with 4 vertices.

julia> minimal_number_of_interior_integral_lines(P)
2
```

"""
function minimal_number_of_interior_integral_lines(P :: RationalPolygon{T}) where {T <: Integer}
    h = ceil(Rational{T}, width(P))
    vs = all_direction_vectors_with_width_less_than(P, h)
    return minimum([number_of_interior_integral_lines(P, v) for v ∈ vs])
end


@doc raw"""
    is_realizable_in_interval(P :: RationalPolygon{T}, h :: T) where {T <: Integer}

Check whether `P` is realizable in ``\mathbb{R} \times [0,h]``. This is true
if and only if `minimal_number_of_interior_integral_lines` is less than or
equal to ``h-1``.

# Example

```jldoctest
julia> P = convex_hull(RationalPoint{Int}[(1//3,-1),(4//3,2),(2//3,2),(-4//3,-1)])
Rational polygon of rationality 3 with 4 vertices.

julia> is_realizable_in_interval(P,2)
false

julia> is_realizable_in_interval(P,3)
true
```

"""
is_realizable_in_interval(P :: RationalPolygon{T}, h :: T) where {T <: Integer} =
minimal_number_of_interior_integral_lines(P) <= h - 1


@doc raw"""
    LatticeWidthData{T <: Integer}

A struct capturing information about the lattice width of a rational polygon
with respect to a some width direction vector, see Definition 2.11 of
[Boh23](@cite). It has two fields:

- `interval_of_nonzero_vertical_slice_length :: Tuple{Rational{T},Rational{T}}`,
- `interval_of_longest_vertical_slice_length :: Tuple{Rational{T},Rational{T}}`.

"""
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


@doc raw"""
    lattice_width_data(P :: RationalPolygon{T}, w :: Point{T}) where {T <: Integer}

Compute the lattice width data of a rational polygon with respect to a given
width direction vector, see Definition 2.11 of [Boh23](@cite). This function returns a value of type `LatticeWidthData`.

# Example

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,1),(1,-2),(-4,2),(-2,2)],2)
Rational polygon of rationality 2 with 4 vertices.

julia> ws = width_direction_vectors(P)
2-element Vector{StaticArraysCore.SVector{2, Int64}}:
 [0, 1]
 [1, 1]

julia> lattice_width_data(P,ws[1])
LatticeWidthData{Int64}((0//1, 2//1), (3//2, 3//2))

julia> lattice_width_data(P,ws[2])
LatticeWidthData{Int64}((0//1, 2//1), (1//2, 1//2))
```

"""
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


@doc raw"""
    number_of_interior_integral_vertical_lines(P :: RationalPolygon, w :: Point{T}) where {T <: Integer}

Return the number of interior integral vertical lines of `P` with respect to a
given lattice width direction vector, see Definition 2.11 of [Boh23](@cite). 

"""
number_of_interior_integral_vertical_lines(P :: RationalPolygon, w :: Point{T}) where {T <: Integer} =
number_of_interior_integral_vertical_lines(lattice_width_data(P, w))


@doc raw"""
    position_of_longest_vertical_slice_length(P :: RationalPolygon, w :: Point{T}) where {T <: Integer}

Return the position of the longest vertical slicing length of `P` with respect
to a given lattice width direction vector, see Definition 2.11 of
[Boh23](@cite). 

"""
position_of_longest_vertical_slice_length(P :: RationalPolygon, w :: Point{T}) where {T <: Integer} =
position_of_longest_vertical_slice_length(lattice_width_data(P, w))


@doc raw"""
    lattice_width_datas(P :: RationalPolygon{T}) where {T <: Integer}   
    
Return the lattice with datas for all lattice width direction vectors of `P`.

"""
lattice_width_datas(P :: RationalPolygon{T}) where {T <: Integer} =
lattice_width_data.(P, width_direction_vectors(P))


@doc raw"""
    numbers_of_interior_integral_vertical_lines(P :: RationalPolygon{T}) where {T <: Integer}

Return the number of interior integral vertical lines for all width direction
vectors of `P`.

"""
numbers_of_interior_integral_vertical_lines(P :: RationalPolygon{T}) where {T <: Integer} =
[number_of_interior_integral_vertical_lines(P, w) for w ∈ width_direction_vectors(P)]


@doc raw"""
    positions_of_longest_vertical_slice_length(P :: RationalPolygon{T}) where {T <: Integer}

Return the positions of the longest vertical slicing lengths for all width
direction vectors of `P`.

"""
positions_of_longest_vertical_slice_length(P :: RationalPolygon{T}) where {T <: Integer} =
[position_of_longest_vertical_slice_length(P, w) for w ∈ width_direction_vectors(P)]
