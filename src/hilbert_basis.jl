
@doc raw"""
    cls_cone_normal_form(A :: Matrix2{T}) where {T <: Integer}

Bring a two-dimensional cone into normal form in the sense of [CLS11](@cite).
The result is a triple (d, k, M), where d and k are the parameters of the cone
and M is a 2x2 integral matrix such that `M * [0 d ; 1 -k] == A`

# Example

```jldoctest
julia> A = Matrix2(2, 1, -3, -5)
2×2 SMatrix{2, 2, Int64, 4} with indices SOneTo(2)×SOneTo(2):
 2  -3
 1  -5

julia> cls_cone_normal_form(A)
(7, 5, [1 2; 0 1])
```

"""
function cls_cone_normal_form(A :: Matrix2{T}) where {T <: Integer}
    (gcd(A[1,1],A[2,1]) == 1 && gcd(A[1,2],A[2,2]) == 1) || error("the given vectors are not primitive")
    # Get integers x.y such that x * A[1,1] + y * A[2,1] == 1 
    _, x, y = gcdx(A[1,1], A[2,1])
    sg, d, l = sign(det(A)), abs(det(A)), x * A[1,2] + y * A[2,2]
    s, k = ((l + mod(-l, d)) ÷ d, mod(-l, d))
    M = @SMatrix [s*A[1,1]-sg*y A[1,1] ; s*A[2,1]+sg*x A[2,1]]
    return (d, k, M)
end


@doc raw"""
    hirzebruch_jung(x :: T, y :: T) where {T <: Integer}

Return the Hirzebruch-Jung continued fraction associated to `x // y`.

# Example

See Example 10.2.4 of [CLS11](@cite).

```jldoctest
julia> hirzebruch_jung(7,5)
3-element Vector{Int64}:
 2
 2
 3
```

"""
function hirzebruch_jung(x :: T, y :: T) where {T <: Integer}
    res = T[]
    while y > 0
        push!(res, (x + mod(-x, y)) ÷ y)
        x, y = y, mod(-x,y)
    end
    return res;
end    


@doc raw"""
    hilbert_basis(A :: Matrix2{T}) where {T <: Integer}

Return the hilbert basis of a two-dimensional cone spanned by the columns of
`A`, which must be primitive.

# Example

See Example 10.2.4 of [CLS11](@cite).

```jldoctest
julia> A = Matrix2(0,1,7,-5)
2×2 SMatrix{2, 2, Int64, 4} with indices SOneTo(2)×SOneTo(2):
 0   7
 1  -5

julia> hilbert_basis(A)
5-element Vector{SVector{2, Int64}}:
 [0, 1]
 [1, 0]
 [2, -1]
 [3, -2]
 [7, -5]
```

"""
function hilbert_basis(A :: Matrix2{T}) where {T <: Integer}
    d, k, M = cls_cone_normal_form(A)

    hj = hirzebruch_jung(d, k)

    x, y = 0, 1
    a, b = -1, 0
    res = Vector{LatticePoint{T}}(undef, length(hj)+2)
    res[1] = M * LatticePoint(0, 1)
    for i = 1 : length(hj)
        z = hj[i]
        res[i+1] = M * LatticePoint(y, -b)
        x, y = y, z * y - x
        a, b = b, z * b - a
    end
    res[end] = M * LatticePoint(y, -b)

    return res
end


@doc raw"""
    hilbert_basis(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}) where {T <: Integer}

Return the hilbert basis of a two-dimensional cone spanned by given integral
primitive vectors `v1` and `v2`.

"""
hilbert_basis(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}) where {T <: Integer} =
hilbert_basis(Matrix2{T}(v1[1],v1[2],v2[1],v2[2]))


function discrepancies(d :: T, k :: T) where {T <: Integer}
    a0, a1 = 1, (k+1) // d
    res = Rational{T}[]
    push!(res, a0)
    for b ∈ hirzebruch_jung(d, k)
        push!(res, a1)
        a0, a1 = a1, b * a1 - a0
    end
    push!(res, a1)
    return res
end

function discrepancies(A :: Matrix2{T}) where {T <: Integer}
    d, k, _ = cls_cone_normal_form(A)
    return discrepancies(d, k)
end


@doc raw"""
    discrepancies(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}) where {T <: Integer}

The discrepancies of the two-dimensional cone spanned by two given integral
primitive vectors.

"""
discrepancies(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}) where {T <: Integer} =
discrepancies(Matrix2{T}(v1[1],v1[2],v2[1],v2[2]))
