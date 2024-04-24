
@attributes mutable struct IntersectionOfHalfplanes{T <: Integer} <: RationalPolygon{T}
    halfplanes :: Vector{AffineHalfplane{T}}

    IntersectionOfHalfplanes(halfplanes :: Vector{AffineHalfplane{T}}) where {T <: Integer} =
    new{T}(halfplanes)

end

Base.:(==)(P1 :: IntersectionOfHalfplanes, P2 :: IntersectionOfHalfplanes) =
P1.halfplanes == P2.halfplanes

affine_halfplanes(P :: IntersectionOfHalfplanes) = P.halfplanes

@attr function vertices(P :: IntersectionOfHalfplanes)
    Hs = affine_halfplanes(P)
    return unique([intersection_point(line(Hs[i]), line(Hs[mod(i+1,1:length(Hs))])) for i = 1 : length(Hs)])
end


