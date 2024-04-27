
@attributes mutable struct IntersectionOfHalfplanes{T <: Integer} <: RationalPolygon{T}
    rationality :: T
    halfplanes :: Vector{AffineHalfplane{T}}

    function IntersectionOfHalfplanes(halfplanes :: Vector{AffineHalfplane{T}}; rationality :: Union{Missing,T} = missing) where {T <: Integer}
        if ismissing(rationality)
            r = length(halfplanes)
            vs = RationalPoint{T}[]
            for i = 1 : r
                H1, H2 = halfplanes[i], halfplanes[mod(i+1,1:r)]
                push!(vs, intersection_point(line(H1),line(H2)))
            end
            unique!(vs)

            k = lcm(RationalPolygons.rationality.(vs))
            P = new{T}(k, halfplanes)
            set_attribute!(P, :vertices, vs)
            return P
        else
            return new{T}(rationality, halfplanes)
        end
    end

end

Base.:(==)(P1 :: IntersectionOfHalfplanes, P2 :: IntersectionOfHalfplanes) =
P1.halfplanes == P2.halfplanes

rationality(P :: IntersectionOfHalfplanes) = P.rationality

affine_halfplanes(P :: IntersectionOfHalfplanes) = P.halfplanes

@attr function vertices(P :: IntersectionOfHalfplanes{T}) where {T <: Integer}
    halfplanes = affine_halfplanes(P)
    r = length(halfplanes)
    vs = RationalPoint{T}[]
    for i = 1 : r
        H1, H2 = halfplanes[i], halfplanes[mod(i+1,1:r)]
        push!(vs, intersection_point(line(H1),line(H2)))
    end
    unique!(vs)
    return vs
end


