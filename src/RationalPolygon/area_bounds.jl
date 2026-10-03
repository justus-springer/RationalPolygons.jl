@doc raw"""
    area_minimizers(k :: T, i :: T, b :: T) where {T <: Integer}

Return the rational polygons of denominator `k ≥ 2` with `i ≥ 1` interior and
`b` boundary lattice points attaining the lower area bounds of Theorem 1.2 in
[BS24_2](@cite), as pairs of their type in that theorem and the polygon. The
polygons of types (0a), (1a), (1b), (2a) and (2b) have the smallest area among
those whose lattice points are not collinear, the polygons of types (0c), (1c)
and (2c) among those whose lattice points are collinear. The returned polygons
are pairwise non-equivalent.

# Example:

```jldoctest
julia> first.(area_minimizers(2, 1, 5))
7-element Vector{String}:
 "1a"
 "2a"
 "2a"
 "2a"
 "2b"
 "2b"
 "2b"
```

"""
function area_minimizers(k :: T, i :: T, b :: T) where {T <: Integer}
    k >= 2 && i >= 1 || throw(ArgumentError("expected k ≥ 2 and i ≥ 1"))
    P(points...) = convex_hull(RationalPoint{T}[points...])
    q = 1 // k
    b_max = (k + 1) * (i + 1) + 3
    Ps = Pair{String,RationalPolygon{T}}[]
    b == b_max - 2(k + 1) && push!(Ps, "0a" => P((0, -1), (i * (k + 1) - k + 1, -1), (-q, q)))
    2 <= b <= b_max - (k + 1) && push!(Ps, "1a" => P((0, -1), (b - 2, -1), (i, 0), (-q, q)))
    (i, b) == (3, 3) && push!(Ps, "1b" => P((0, -2), (2, 0), (-q, q)))
    3 <= b <= b_max && append!(Ps, ["2a" => P((0, 0), (0, -1), (b - 3, -1), (i + 1, 0), (x * q, q))
                                    for x in 0 : (b_max - b) ÷ 2])
    (i, b) == (1, 5) && append!(Ps, ["2b" => P((0, 0), (0, -2), (2, 0), (x * q, q)) for x in 0 : k])
    b == 0 && push!(Ps, "0c" => P((1, q), (1 - q, -q), (i + q, 0)))
    b == 1 && push!(Ps, "1c" => P((1, q), (1 - q, -q), (i + 1, 0)))
    b == 2 && append!(Ps, ["2c" => P((0, 0), (0, q), (x * q, -q), (i + 1, 0)) for x in 0 : k * (i + 1)])
    return Ps
end


@doc raw"""
    area_maximizers(k :: T, i :: T, b :: T) where {T <: Integer}

Return the rational polygons of denominator `k ≥ 2` with `i ≥ 1` interior and
`b` boundary lattice points attaining the upper area bound of Theorem 1.3 in
[BS24_2](@cite), as pairs of their type in that theorem and the polygon. For
`k ≥ 4`, these are precisely the polygons of maximal area. For `k = 3`, this
has been verified for `i ≤ 5`, see Remark 2.7 in [BS24_2](@cite). For `k = 2`,
the bound differs if `i = 1` or `b = 0` and there are further maximizers, see
[`maximal_area_half_integral`](@ref). The returned polygons are pairwise
non-equivalent.

# Example:

```jldoctest
julia> first.(area_maximizers(4, 1, 2))
6-element Vector{String}:
 "1a"
 "2a"
 "2a"
 "2a"
 "2a"
 "2a"
```

"""
function area_maximizers(k :: T, i :: T, b :: T) where {T <: Integer}
    k >= 2 && i >= 1 || throw(ArgumentError("expected k ≥ 2 and i ≥ 1"))
    P(points...) = convex_hull(RationalPoint{T}[points...])
    q = 1 // k
    b_max = (k + 1) * (i + 1) + 3
    Ps = Pair{String,RationalPolygon{T}}[]
    b == b_max - 4 && push!(Ps, "0a" => P((0, q), (q, -1), ((k + 1) * (i + 1) - q, -1)))
    b == 0 && push!(Ps, "0b" => P((0, q), (q, -1), (1 - q, -1), (k * (i + 1) - q, q - 1)))
    1 <= b <= b_max - 3 && push!(Ps, "1a" => P((0, q), (q, -1), (b - q, -1), (k * (i + 1), q - 1)))
    2 <= b <= b_max - 2 && append!(Ps, ["2a" => P((0, q), (0, q - 1), (x + q, -1), (x + b - 1 - q, -1), (k * (i + 1), q - 1))
                                        for x in 0 : (b_max - b) ÷ 2 - 1])
    b == b_max - 1 && push!(Ps, "2b" => P((0, q), (0, q - 1), (q, -1), ((k + 1) * (i + 1), -1)))
    b == b_max && push!(Ps, "2c" => P((0, q), (0, -1), ((k + 1) * (i + 1), -1)))
    return Ps
end


@doc raw"""
    maximal_area_half_integral(i :: T, b :: T) where {T <: Integer}

Return the maximal euclidean area of a rational polygon of denominator two
with `i ≥ 1` interior and `b` boundary lattice points, see Theorem 1.4 in
[BS24_2](@cite).

# Example:

```jldoctest
julia> maximal_area_half_integral(1, 0)
21//8
```

"""
function maximal_area_half_integral(i :: T, b :: T) where {T <: Integer}
    i >= 1 && 0 <= b <= 3i + 6 || throw(ArgumentError("expected i ≥ 1 and 0 ≤ b ≤ 3i + 6"))
    i == 1 && return b // 4 + (b <= 6 ? 21 : 27 - b) // 8
    return 3i // 2 + b // 4 + (b <= 3i + 4 ? 8 : b == 3i + 5 ? 7 : 6) // 8
end
