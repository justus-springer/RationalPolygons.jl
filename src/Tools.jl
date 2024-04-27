
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
