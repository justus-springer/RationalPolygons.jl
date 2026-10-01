@doc raw"""
    hnf(A :: SMatrix{2,N,T}) where {N,T<:Integer}

An elementary implementation of the hermite normal form for 2xn integral
matrices.

"""
function hnf(A :: SMatrix{2,N,T,M}) where {N,M,T<:Integer}

    # The pivots are the first nonzero column and the first column independent
    # of it. These need not be the first two columns, e.g. for the vertex
    # matrix of a polygon with the origin on an edge line.
    i = 1
    while A[1,i] == 0 && A[2,i] == 0
        i == N && return A
        i += 1
    end

    if A[1,i] == 0
        A = SMatrix{2,2,T,4}(0,1,1,0) * A
    end

    d,a,b = gcdx(A[1,i],A[2,i])
    if A[1,i] != d
        _, x, y = gcdx(a,b)
        A = SMatrix{2,2,T,4}(a,-y,b,x) * A
    end

	f = div(A[2,i],A[1,i])
    A = SMatrix{2,2,T,4}(1,-f,0,1) * A

    j = i + 1
    while j <= N && A[2,j] == 0
        j += 1
    end
    j > N && return A

    if sign(A[2,j]) == -1
        A = SMatrix{2,2,T,4}(1,0,0,-1) * A
    end

	c = fld(A[1,j],A[2,j])
    if c != 0
        A = SMatrix{2,2,T,4}(1,0,-c,1) * A
    end

    return A

end


@doc raw"""
    lattice_edge_areas(P :: RationalPolygon)

To each vertex of a polygon, we associate the area spanned by the two adjacent
edges. This function returns the vector of all those areas.

"""
lattice_edge_areas(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
SVector{N}(abs(det(P[i+1] - P[i], P[i] - P[i-1])) for i = 1 : N)


@doc raw"""
    area_maximizing_vertices(P :: RationalPolygon)

Return the indices of those vertices of `P` that maximize the lattice edge
area.

"""
function area_maximizing_vertices(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    ea = lattice_edge_areas(P)
    m = maximum(ea)
    return filter(i -> ea[i] == m, 1 : N)
end

function _translate_columns(A :: SMatrix{2,N,T,M}, i :: Int, reverse_columns :: Bool) where {N,M,T}
    if reverse_columns
        return SMatrix{2,N,T,2N}(isodd(j) ? A[1,mod(i-(j+1)÷2,1:N)] : A[2,mod(i-(j+1)÷2,1:N)] for j = 1 : 2N)
    else
        return SMatrix{2,N,T,2N}(isodd(j) ? A[1,mod(i+(j+1)÷2,1:N)] : A[2,mod(i+(j+1)÷2,1:N)] for j = 1 : 2N)
    end

end

@doc raw"""
    unimodular_normal_form_with_automorphism_group(P :: RationalPolygon)

Return a pair `(Q,G)` where `Q` is the unimodular normal form of `P` and `G` is
the unimodular automorphism group of `P`.

"""
function unimodular_normal_form_with_automorphism_group(P :: RationalPolygon{T,N}) where {N,T <: Integer}

    V = vertex_matrix(P)

    # For all area maximizing vertices `v`, we consider the vertex numbering
    # starting at `v` and going clockwise or counterclockwise. For all these
    # possible numberings, we compute the hermite normal form of the vertex 
    # matrix
    
    As = SVector{2N,SMatrix{2,N,T,2N}}(hnf(_translate_columns(V,(i÷2)+1,isodd(i))) for i = 0 : 2N-1)
    A = argmin(vec, As)

    special_indices = map(l -> divrem(l+1, 2), findall(B -> A == B, As))

    Q = RationalPolygon(A, rationality(P); is_unimodular_normal_form = true) 
    if all(sv -> sv[2] == 0, special_indices) || all(sv -> sv[2] == 1, special_indices)
        return (Q, CyclicGroup(length(special_indices)))
    else
        return (Q, DihedralGroup(length(special_indices) ÷ 2))
    end

end

@doc raw"""
    unimodular_normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer}

Return a unimodular normal form of a rational polygon. Two rational polygons have
the same unimodular normal form if and only if the can be transformed into each
other by applying a unimodular transformation.

"""
function unimodular_normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    is_unimodular_normal_form(P) && return P

    V = vertex_matrix(P)

    As = SVector{2N,SMatrix{2,N,T,2N}}(hnf(_translate_columns(V,(i÷2)+1,isodd(i))) for i = 0 : 2N-1)

    A = argmin(vec, As)
    Q = RationalPolygon(A, rationality(P); is_unimodular_normal_form = true) 

    return Q

end


@doc raw"""
    prettify(P :: RationalPolygon{T}) where {T <: Integer}

Return a unimodular image of `P` whose bounding box has minimal side lengths.

This uses the generalized Gauss reduction, see Algorithm 2.5 of
[HaSo22](@cite): a reduced basis with respect to the width norm
`h ↦ width(P, h)` attains both successive minima, see Theorem 2.3 of
[HaSo22](@cite).

"""
function prettify(P :: RationalPolygon{T}) where {T <: Integer}
    V = vertex_matrix(P)
    nrm(h) = maximum(h' * V) - minimum(h' * V)

    h1, h2 = SVector{2,T}(1,0), SVector{2,T}(0,1)
    # Ties are broken towards the identity, so that prettify is idempotent.
    nrm(h1) >= nrm(h2) && ((h1, h2) = (h2, h1))
    while true
        # m ↦ nrm(m*h1 + h2) is convex and its minimum lies in [-M, M]
        # by the triangle inequality, so binary search for the first
        # non-negative difference.
        M = cld(2 * nrm(h2), nrm(h1))
        lo, hi = -M, M
        while lo < hi
            m = fld(lo + hi, 2)
            nrm((m+1) * h1 + h2) >= nrm(m * h1 + h2) ? (hi = m) : (lo = m + 1)
        end
        h = nrm(h2) == nrm(lo * h1 + h2) ? h2 : lo * h1 + h2
        if nrm(h) >= nrm(h1)
            h2 = h
            break
        end
        h1, h2 = h, h1
    end

    return vcat(h2', h1') * P
end


@doc raw"""
    pretty_normal_form(P :: RationalPolygon{T}) where {T <: Integer}

Return `prettify(unimodular_normal_form(P))`, a canonical representative with a
minimal bounding box.

"""
pretty_normal_form(P :: RationalPolygon{T}) where {T <: Integer} =
prettify(unimodular_normal_form(P))


@doc raw"""
    are_unimodular_equivalent(P :: RationalPolygon, Q :: RationalPolygon)   

Checks whether two rational polygons are equivalent by a unimodular
transformation.

"""
are_unimodular_equivalent(P :: RationalPolygon, Q :: RationalPolygon) =
unimodular_normal_form(P) == unimodular_normal_form(Q)

@doc raw"""
    unimodular_automorphism_group(P :: RationalPolygon)

Return the automorphism group of `P` with respect to unimodular
transformations.

# Example:

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)])
Rational polygon of rationality 1 with 6 vertices.

julia> unimodular_automorphism_group(P)
D6
```

"""
unimodular_automorphism_group(P :: RationalPolygon) =
unimodular_normal_form_with_automorphism_group(P)[2]


@doc raw"""
    function _special_vertices_align(k :: T,
        A1 :: SMatrix{2,N,T}, i1 :: Int, o1 :: Bool,
        A2 :: SMatrix{2,N,T}, i2 :: Int, o2 :: Bool) where {N, T <: Integer}

This is a helper function for `affine_normal_form`. It checks, whether two
special vertices of two `k`-rational polygons can be sent to each other via an
affine unimodular transformation. The matrices `A1, A2` are the vertex matrices
of the polygons, `i1, i2` are the indices of the special vertices and `o1, o2`
are the orientations, i.e. whether the vertices of `A1` and `A2` are sorted
clockwise or counterclockwise.

"""
function _special_vertices_align(k :: T,
        A1 :: SMatrix{2,N,T}, i1 :: Int, o1 :: Int,
        A2 :: SMatrix{2,N,T}, i2 :: Int, o2 :: Int) where {N, T <: Integer}
    s1 = o1 == 0 ? 1 : -1
    s2 = o2 == 0 ? 1 : -1
    v1, v2, v3 = A1[:,mod(i1-s1,1:N)], A1[:,mod(i1,1:N)], A1[:,mod(i1+s1,1:N)]
    w1, w2, w3 = A2[:,mod(i2-s2,1:N)], A2[:,mod(i2,1:N)], A2[:,mod(i2+s2,1:N)]
    d1, d2, d3 = det(w2,w3), det(w3,w1), det(w1,w2)
    d = d1 + d2 + d3
    b1 = (v1[1] * d1 + v2[1] * d2 + v3[1] * d3) // d
    b2 = (v1[2] * d1 + v2[2] * d2 + v3[2] * d3) // d
    return b1 % k == 0 && b2 % k == 0
end

@doc raw"""
    affine_normal_form_with_special_indices(P :: RationalPolygon)

Return a pair `(Q,is)` where `Q` is the affine normal form of `P` and `is` is the set of special indices of `P`.

"""
function affine_normal_form_with_special_indices(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    V = vertex_matrix(P)
    k = rationality(P)

    As = SVector{2N,SMatrix{2,N,T,2N}}(hnf(_translate_columns(V,(i÷2)+1,isodd(i)).- scaled_vertex(P,(i÷2)+1)) for i = 0 : 2N-1)
    A = argmin(vec, As)

    special_indices = Tuple{Int,Int}[]
    for l = 1 : 2N
        if As[l] == A
            push!(special_indices, divrem(l+1, 2))
        end
    end

    really_special_indices = Tuple{Int,Int}[]
    local At :: SMatrix{2,N,T,2N}
    for x = 0 : k-1, y = 0 : k-1
        At = A .+ LatticePoint{T}(x,y)
        empty!(really_special_indices)
        for (j,o) ∈ special_indices
            if _special_vertices_align(k, At, N, 0, V, j, o)
                push!(really_special_indices, (j,o))
            end
        end
        !isempty(really_special_indices) && break
    end

    Q = RationalPolygon(At, rationality(P); is_affine_normal_form = true)

    return (Q, really_special_indices)

end


@doc raw"""
    affine_normal_form_with_automorphism_group(P :: RationalPolygon)

Return a pair `(Q,G)` where `Q` is the affine normal form of `P` and `G` is the
affine automorphism group of `P`.

"""
function affine_normal_form_with_automorphism_group(P :: RationalPolygon{T,N}) where {N,T <: Integer}

    (Q,is) = affine_normal_form_with_special_indices(P)

    if length(unique(last.(is))) == 1
        return (Q, CyclicGroup(length(is)))
    else
        return (Q, DihedralGroup(length(is) ÷ 2))
    end

end


@doc raw"""
    affine_normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer}

Return an affine normal form of a rational polygon. Two rational polygons have
the same affine normal form if and only if the can be transformed into each
other by applying an affine unimodular transformation.

"""
affine_normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
is_affine_normal_form(P) ? P : affine_normal_form_with_automorphism_group(P)[1]


@doc raw"""
    are_affine_equivalent(P :: RationalPolygon, Q :: RationalPolygon)   

Checks whether two rational polygons are equivalent by an affine unimodular
transformation.

"""
are_affine_equivalent(P :: RationalPolygon, Q :: RationalPolygon) =
affine_normal_form(P) == affine_normal_form(Q)

@doc raw"""
    affine_automorphism_group(P :: RationalPolygon{T,N}) where {N,T <: Integer}

Return the automorphism group of `P` with respect to affine unimodular
transformations.

"""
affine_automorphism_group(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
affine_normal_form_with_automorphism_group(P)[2]


function lattice_affine_normal_form_with_is_first_index_special(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    V = vertex_matrix(P)

    As = SVector{2N,SMatrix{2,N,T,2N}}(hnf(_translate_columns(V,(i÷2)+1,isodd(i)).- scaled_vertex(P,(i÷2)+1)) for i = 0 : 2N-1)
    A = argmin(vec, As)
    Q = RationalPolygon(A, one(T); is_affine_normal_form = true)

    is_first_index_special = As[1] == A || As[2] == A

    return (Q, is_first_index_special)

end
