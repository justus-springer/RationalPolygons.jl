
floor_k_rational(k :: T, x :: Real) where {T <: Integer} =
floor(T, k * x) // k

ceil_k_rational(k :: T, x :: Real) where {T <: Integer} =
ceil(T, k * x) // k

function pseudo_angle(p :: Point{T}) where {T <: Integer}
    x,y = p[1],p[2]
    (x == 0 && y == 0) && return 0
    d = abs(x) + abs(y)
    a = x // d 
    y < 0 && return a - 1
    return 1 - a
end

function pseudo_angle_with_distance(p :: Point{T}) where {T <: Integer}
    x,y = p[1],p[2]
    (x == 0 && y == 0) && return (0,0)
    d = abs(x) + abs(y)
    a = x // d 
    y < 0 && return (a - 1, d)
    return (1 - a, d)
end

# make Julia use det_bareiss by default for integral matrices
det(M :: AbstractMatrix{T}) where {T <: Integer} = det_bareiss(M)

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


