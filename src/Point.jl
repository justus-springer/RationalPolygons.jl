
const LatticePoint{T<:Integer} = Tuple{T,T}

const RationalPoint{T<:Integer} = Tuple{Rational{T},Rational{T}}

const Point{T<:Integer} = Tuple{S,S} where S<:Union{T,Rational{T}}

rationality(p :: Point) = lcm(denominator(p[1]), denominator(p[2]))
multiplicity(p :: LatticePoint) = gcd(p[1], p[2])

numerator(p :: Point{T}) where {T <: Integer} = numerator.(p)
denominator(p :: Point{T}) where {T <: Integer} = denominator.(p)
Base.:(+)(p :: Point{T}, q :: Point{T}) where {T <: Integer} = p .+ q
Base.:(*)(c :: Union{T,Rational{T}}, p :: Point{T}) where {T <: Integer} = c .* p
Base.:(//)(p :: Point{T}, c :: Union{T,Rational{T}}) where {T <: Integer} = p .// c
