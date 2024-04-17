
struct RationalPolygon{T<:Integer}
    rationality :: T
    vs :: Vector{Tuple{T,T}}

    RationalPolygon(rationality :: T, vs :: Vector{Tuple{T,T}}) where {T <: Integer} =
    new{T}(rationality, vs)

    function RationalPolygon(vs :: Vector{Tuple{T,T}}) where {T <: Union{Integer,Rational}}
        r = lcm([lcm(denominator(v[1]), denominator(v[2])) for v ∈ vs])
        ws = [(r * numerator(v[1]) ÷ denominator(v[1]), r * numerator(v[2]) ÷ denominator(v[2])) for v ∈ vs]
        return new{typeof(r)}(r,ws)
    end

end

Base.show(io :: IO, P :: RationalPolygon) =
Base.print(io, "Rational polygon of rationality $(P.rationality)")

Base.:(==)(P :: RationalPolygon, Q :: RationalPolygon) =
P.rationality == Q.rationality && P.vs == Q.vs

rationality(P :: RationalPolygon) = P.rationality

scaled_vertices(P :: RationalPolygon) = P.vs

vertices(P :: RationalPolygon) = [v .// rationality(P) for v ∈ scaled_vertices(P)]
