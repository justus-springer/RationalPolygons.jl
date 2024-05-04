
@doc raw"""
    cls_cone_normal_form(A :: Matrix{T}) where {T <: Integer}

Bring a two-dimensional cone into normal form in the sense of Cox, Little,
Schenck Proposition 10.1.1. The result is a triple (d, k, M), where d and k are
the parameters of the cone and M is a 2x2 integral matrix such that 
`M * [0 d ; 1 -k] == A`

"""
function cls_cone_normal_form(A :: Matrix{T}) where {T <: Integer}
    u11, u21, u12, u22 = A[1,1], A[2,1], A[1,2], A[2,2]
    (gcd(u11,u21) == 1 && gcd(u12,u22) == 1) || error("the given vectors are not primitive")
    # Get integers x.y such that x * u11 + y * u21 == 1 
    _, x, y = gcdx(u11, u21)
    sg, d, l = sign(det(A)), abs(det(A)), x * u12 + y * u22
    s, k = ((l + mod(-l, d)) ÷ d, mod(-l, d))
    M = [s*u11-sg*y u11 ; s*u21+sg*x u21]
    return (d, k, M)
end


@doc raw"""
    hirzebruch_jung(x :: T, y :: T)

Return the Hirzebruch-Jung continued fraction associated to `x // y`.

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
    hilbert_basis(A :: Matrix{T}) where {T <: Integer}

Return the hilbert basis of a two-dimensional cone spanned by the columns of
`A`, which must be primitive.

"""
function hilbert_basis(A :: Matrix{T}) where {T <: Integer}
    d, k, M = cls_cone_normal_form(A)

    hj = hirzebruch_jung(d, k)

    x, y = 0, 1
    a, b = -1, 0
    res = Vector{T}[]
    push!(res, M * [0 ; 1])
    for z in hj
        push!(res, M * [y ; -b])
        x, y = y, z * y - x
        a, b = b, z * b - a
    end
    push!(res, M * [y ; -b])

    res = map(v -> (v[1],v[2]), res)
    return res
end


@doc raw"""
    hilbert_basis(v1 :: Vector{T}, v2 :: Vector{T})

Return the hilbert basis of a two-dimensional cone spanned by given integral
primitive vectors `v1` and `v2`.

"""
hilbert_basis(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}) where {T <: Integer} =
hilbert_basis([v1[1] v2[1] ; v1[2] v2[2]])


