
# an elementary implementation of the hermite normal form for 2xn matrices
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

function special_vertices(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    ea = lattice_edge_areas(P)
    m = maximum(ea)
    return filter(i -> ea[i] == m, 1 : N)
end

function normal_form(P :: RationalPolygon{T,N}) where {N,T <: Integer}

    is_normal_form(P) && return P

    V = vertex_matrix(P)

    As = MMatrix{2,N,T,2N}[]
    for i ∈ special_vertices(P)
        push!(As, MMatrix{2,N,T,2N}([V[:,i:end] V[:,begin:i-1]]))
        push!(As, MMatrix{2,N,T,2N}([V[:,i:-1:begin] V[:,end:-1:i+1]]))
    end
    hnf!.(As)

    # take the lexicographical minumum of all the hermite normal forms
    _lt(A, B) = vec(A) < vec(B)
    A = sort(As; lt = _lt)[1]

    Q = RationalPolygon(convert(SMatrix,A), rationality(P); is_normal_form = true) 

    return Q

end

function are_equivalent(P :: RationalPolygon, Q :: RationalPolygon)
    return normal_form(P) == normal_form(Q)
end
