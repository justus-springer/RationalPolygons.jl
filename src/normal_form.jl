
@doc raw"""
    hnf!(A :: MMatrix{2,N,T}) where {N,T<:Integer}

An elementary implementation of the hermite normal form for 2xn integral
matrices.

"""
function hnf!(A :: MMatrix{2,N,T}) where {N,T<:Integer}

    A[1,1] == 0 && mul!(A, SMatrix{2,2,T,4}(0,1,1,0), copy(A))

    d,a,b = gcdx(A[1,1],A[2,1])
    if A[1,1] != d
        _, x, y = gcdx(a,b)
        mul!(A, SMatrix{2,2,T,4}(a,-y,b,x), copy(A))
    end

	if A[1,1] != 0 
		f = div(A[2,1],A[1,1])
        mul!(A, SMatrix{2,2,T,4}(1,-f,0,1), copy(A))
	end

    sign(A[2,2]) == -1 && mul!(A, SMatrix{2,2,T,4}(1,0,0,-1), copy(A))

	if A[2,2] != 0
		c = fld(A[1,2],A[2,2])
        c != 0 && mul!(A, SMatrix{2,2,T,4}(1,0,-c,1), copy(A))
	end

    return A

end

lattice_edge_areas(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
[abs(det(P[i+1] - P[i], P[i] - P[i-1])) for i = 1 : N]

function area_maximizing_vertices(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    ea = lattice_edge_areas(P)
    m = maximum(ea)
    return filter(i -> ea[i] == m, 1 : N)
end

function unimodular_normal_form_with_automorphism_group(P :: RationalPolygon{T,N}) where {N,T <: Integer}

    V = vertex_matrix(P)

    # For all area maximizing vertices `v`, we consider the vertex numbering
    # starting at `v` and going clockwise or counterclockwise. For all these
    # possible numberings, we compute the hermite normal form of the vertex 
    # matrix
    As = MMatrix{2,N,T,2N}[]
    for i ∈ area_maximizing_vertices(P)
        push!(As, MMatrix{2,N,T,2N}([V[:,i:end] V[:,begin:i-1]]))
        push!(As, MMatrix{2,N,T,2N}([V[:,i:-1:begin] V[:,end:-1:i+1]]))
    end
    hnf!.(As)

    # We take the normal form to be the lexicographical minimum of all those
    # hermite normal forms
    A = argmin(vec, As)
    special_indices = map(l -> divrem(l+1, 2), findall(B -> A == B, As))

    Q = RationalPolygon(convert(SMatrix,A), rationality(P); is_unimodular_normal_form = true) 
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
unimodular_normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
is_unimodular_normal_form(P) ? P : unimodular_normal_form_with_automorphism_group(P)[1]


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


function affine_normal_form_with_automorphism_group(P :: RationalPolygon{T,N}) where {N,T <: Integer}

    V = vertex_matrix(P)
    k = rationality(P)

    As = MMatrix{2,N,T,2N}[]
    is = area_maximizing_vertices(P)

    for i ∈ is
        # move the vertex to the origin
        push!(As, MMatrix{2,N,T,2N}([V[:,i+1:end] V[:,begin:i]]) .- lattice_vertex(P,i))
        push!(As, MMatrix{2,N,T,2N}([V[:,i-1:-1:begin] V[:,end:-1:i]]) .- lattice_vertex(P,i))
    end
    hnf!.(As)

    A = argmin(vec, As)
    special_indices = Tuple{Int,Int}[]
    for l = 1 : length(As)
        if As[l] == A
            push!(special_indices, divrem(l+1, 2))
        end
    end
    A = convert(SMatrix,A)

    really_special_indices = Tuple{Int,Int}[]
    local At
    for x = 0 : k-1, y = 0 : k-1
        At = A .+ LatticePoint{T}(x,y)
        empty!(really_special_indices)
        for (j,o) ∈ special_indices
            if _special_vertices_align(k, At, N, 0, V, is[j], o)
                push!(really_special_indices, (j,o))
            end
        end
        !isempty(really_special_indices) && break
    end

    Q = RationalPolygon(At, rationality(P); is_affine_normal_form = true)

    if length(unique(last.(really_special_indices))) == 1
        return (Q, CyclicGroup(length(really_special_indices)))
    else
        return (Q, DihedralGroup(length(really_special_indices) ÷ 2))
    end

end


@doc raw"""
    affine_normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer}

Return a affine normal form of a rational polygon. Two rational polygons have
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
