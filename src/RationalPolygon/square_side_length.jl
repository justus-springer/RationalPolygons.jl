@doc raw"""
    square_side_length(P :: RationalPolygon{T}) where {T <: Integer}

Return the smallest integer ``m`` such that `P` is affine unimodular
equivalent to a subpolygon of the square ``[0,m]^2``. For lattice polygons,
these are classified by Brown and Kasprzyk [BK13](@cite) and in Table 2.13 of
[BS24_1](@cite).

# Example:

```jldoctest
julia> square_side_length(convex_hull(LatticePoint{Int}[(0,0),(3,0),(0,3)]))
3
```

"""
function square_side_length(P :: RationalPolygon{T}) where {T <: Integer}
    V, k = vertex_matrix(P), rationality(P)
    # P lies between the integral lines <w,x> = c and <w,x> = c + m for m at least this.
    side(w) = cld(maximum(w[1] * V[1,j] + w[2] * V[2,j] for j in axes(V, 2)), k) -
              fld(minimum(w[1] * V[1,j] + w[2] * V[2,j] for j in axes(V, 2)), k)
    # Double c, starting from the lower bound ⌈√area⌉, until two directions with side at most c form a basis. The
    # optimal basis then consists of directions of width at most c. Searching only directions of small width keeps this
    # fast, whereas the standard basis or the lattice width can be costly for skewed polygons.
    c = max(one(T), T(isqrt(cld(numerator(euclidean_area(P)), denominator(euclidean_area(P))))))
    while true
        ws = filter(w -> side(w) <= c, all_direction_vectors_with_width_less_than(P, c // one(T)))
        sides = side.(ws)
        m = minimum((max(sides[a], sides[b]) for a in eachindex(ws), b in eachindex(ws)
                     if abs(ws[a][1] * ws[b][2] - ws[a][2] * ws[b][1]) == 1); init = c + 1)
        m <= c && return m
        c *= 2
    end
end
