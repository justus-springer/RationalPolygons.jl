
const LatticePoint{T<:Integer} = Tuple{T,T}

const RationalPoint{T<:Integer} = Tuple{Rational{T},Rational{T}}

const Point{T<:Integer} = Tuple{S,S} where S<:Union{T,Rational{T}}

rationality(p :: Point) = lcm(denominator(p[1]), denominator(p[2]))
multiplicity(p :: LatticePoint) = gcd(p[1], p[2])

numerator(p :: Point{T}) where {T <: Integer} = numerator.(p)
denominator(p :: Point{T}) where {T <: Integer} = denominator.(p)
Base.:(+)(p :: Point{T}, q :: Point{T}) where {T <: Integer} = p .+ q
Base.:(-)(p :: Point{T}, q :: Point{T}) where {T <: Integer} = p .- q
Base.:(-)(p :: Point{T}) where {T <: Integer} = .-p
Base.:(*)(c :: Union{T,Rational{T}}, p :: Point{T}) where {T <: Integer} = c .* p
Base.:(//)(p :: Point{T}, c :: Union{T,Rational{T}}) where {T <: Integer} = p .// c
Base.zero(::LatticePoint{T}) where {T <: Integer} = (zero(T), zero(T))
Base.zero(::Type{LatticePoint{T}}) where {T <: Integer} = (zero(T), zero(T))
Base.zero(::RationalPoint{T}) where {T <: Integer} = (zero(T) // one(T), zero(T) // one(T))
Base.zero(::Type{RationalPoint{T}}) where {T <: Integer} = (zero(T) // one(T), zero(T) // one(T))

norm(p :: Point{T}) where {T <: Integer} = p[1]^2 + p[2]^2
distance(p :: Point{T}, q :: Point{T}) where {T <: Integer} = norm(p - q)

det(p :: Point{T}, q :: Point{T}) where {T <: Integer} = p[1] * q[2] - p[2] * q[1]

is_k_rational(k :: T, p :: Point{T}) where {T <: Integer} = k % rationality(p) == 0

is_integral(p :: Point{T}) where {T <: Integer} = is_k_rational(one(T), p)
