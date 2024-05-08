
const LatticePoint{T <: Integer} = SVector{2, T}

const RationalPoint{T<:Integer} = SVector{2, Rational{T}}

const Point{T<:Integer} = SVector{2, S} where S<:Union{T,Rational{T}}

const Matrix2{T <: Integer} = SMatrix{2, 2, T, 4}

det(p :: Point{T}, q :: Point{T}) where {T <: Integer} = p[1] * q[2] - p[2] * q[1]

det(M :: Matrix2{T}) where {T <: Integer} = M[1,1]*M[2,2] - M[1,2]*M[2,1]


@doc raw"""
    rationality(p :: Point)

The smallest integer `r` such that `r*p` is integral.

"""
rationality(p :: Point) = lcm(denominator(p[1]), denominator(p[2]))


@doc raw"""
    multiplicity(p :: Point)

The unique rational number `x` such that `x*p` is primitive and integral.

"""
multiplicity(p :: Point) = (denominator(p[1]) * denominator(p[2])) // gcd(numerator(p[1]) * denominator(p[2]), numerator(p[2]) * denominator(p[1]))

is_primitive(p :: Point) = multiplicity(p) == 1
primitivize(p :: Point) = numerator.(multiplicity(p) * p)

norm(p :: Point{T}) where {T <: Integer} = p[1]^2 + p[2]^2
distance(p :: Point{T}, q :: Point{T}) where {T <: Integer} = norm(p - q)

is_k_rational(k :: T, p :: Point{T}) where {T <: Integer} = k % rationality(p) == 0

is_integral(p :: Point{T}) where {T <: Integer} = is_k_rational(one(T), p)

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
