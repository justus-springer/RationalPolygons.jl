
@doc raw"""
    LatticePoint{T <: Integer}

A lattice point in two-dimensional space. This is an alias for `SVector{2, T}`.

"""
const LatticePoint{T <: Integer} = SVector{2, T}


@doc raw"""
    RationalPoint{T<:Integer}

A rational point in two-dimensional space. This is an alias for `SVector{2,
Rational{T}}`.

"""
const RationalPoint{T<:Integer} = SVector{2, Rational{T}}


@doc raw"""
    Point{T<:Integer} 

The union of `LatticePoint` and `RationalPoint`.

"""
const Point{T<:Integer} = SVector{2, S} where S<:Union{T,Rational{T}}


@doc raw"""
    Matrix2{T <: Integer}

The type of 2x2 matrices over the integers 

"""
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


@doc raw"""
    is_primitive(p :: Point)

Checks whether a point is primitive, i.e. is integral and its coordinates are
coprime.

"""
is_primitive(p :: Point) = multiplicity(p) == 1


@doc raw"""
    primitivize(p :: Point)

Return the unique primitive lattice point on the ray spanned by `p`.

"""
primitivize(p :: Point) = numerator.(multiplicity(p) * p)

numerator(p :: Point) = numerator.(p)
denominator(p :: Point) = denominator.(p)


@doc raw"""
    norm(p :: Point{T})

Return the square of the euclidian norm of `p`.

"""
norm(p :: Point{T}) where {T <: Integer} = p[1]^2 + p[2]^2


@doc raw"""
    distance(p :: Point{T}, q :: Point{T})

Return the square of the euclidian distance between `p` and `q`.

"""
distance(p :: Point{T}, q :: Point{T}) where {T <: Integer} = norm(p - q)


@doc raw"""
    is_k_rational(k :: T, p :: Point{T}) where {T <: Integer}

Checks whether a point `p` is `k`-rational, i.e. its coordinates have
denominator at most `k`.

"""
is_k_rational(k :: T, p :: Point{T}) where {T <: Integer} = k % rationality(p) == 0


@doc raw"""
    is_integral(p :: Point{T}) where {T <: Integer}

Checks whether a point `p` is integral.

"""
is_integral(p :: Point{T}) where {T <: Integer} = is_k_rational(one(T), p)


@doc raw"""
    pseudo_angle(p :: Point{T}) where {T <: Integer}

A quick implementation of a pseudo_angle of two-dimensional vectors. It
returns values in the half-open interval (-2,2].

"""
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

norm_ratio(p :: Point{T}, q :: Point{T}) where {T <: Integer} =
q[2] == 0 ? p[1] // q[1] : p[2] // q[2]
