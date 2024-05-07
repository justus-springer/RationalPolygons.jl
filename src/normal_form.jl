
# an elementary implementation of the hermite normal form for 2xn matrices
function hnf!(A :: AbstractMatrix{<:Integer})
    n = size(A,2)

    A[1,1] == 0 && (A = [0 1 ; 1 0] * A)
    d,a,b = gcdx(A[1,1],A[2,1])

    if A[1,1] != d
        _, x, y = gcdx(a,b)
        A = [a b ; -y x] * A
    end

	if A[1,1] != 0 
		f = div(A[2,1],A[1,1])
        A = [1 0 ; -f 1] * A
	end

    sign(A[2,2]) == -1 && (A = [1 0 ; 0 -1] * A)

	if A[2,2] != 0
		c = fld(A[1,2],A[2,2])
        c != 0 && (A = [1 -c ; 0 1] * A)
	end

    return A

end

hnf(A :: AbstractMatrix{<:Integer}) = hnf!(copy(A))


function lattice_edge_areas(P :: RationalPolygon)
    n = number_of_vertices(P)
    return [abs(det(P[i+1] - P[i], P[i] - P[i-1])) for i = 1 : n]
end

function special_vertices(P :: RationalPolygon)
    r = number_of_vertices(P)
    ea = lattice_edge_areas(P)
    m = maximum(ea)
    return filter(i -> ea[i] == m, 1 : r)
end

function normal_form(P :: RationalPolygon{T}) where {T <: Integer}

    V = vertex_matrix(P)

    As = Matrix{T}[]
    for i ∈ special_vertices(P)
        push!(As, [V[:,i:end] V[:,begin:i-1]])
        push!(As, [V[:,i:-1:begin] V[:,end:-1:i+1]])
    end
    As = map(hnf, As)

    # take the lexicographical minumum of all the hermite normal forms
    _lt(A, B) = vec(A) < vec(B)
    A = sort(As; lt = _lt)[1]

    Q = RationalPolygon(A, rationality(P)) 

    return Q

end

is_normal_form(P :: RationalPolygon) =
P == normal_form(P)

function are_equivalent(P :: RationalPolygon, Q :: RationalPolygon)
    return normal_form(P) == normal_form(Q)
end
