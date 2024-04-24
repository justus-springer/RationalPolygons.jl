
@doc raw"""
    abstract type AffineHalfplane end

An affine halfplane in two-dimensional rational space. Subtypes must at least
implement `normal_vector` and `translation`.

"""
abstract type AffineHalfplane{T <: Integer} end

function Base.show(io :: IO, H :: AffineHalfplane)
    nv, b = normal_vector(H), translation(H)
    print(io, "Affine halfplane given by $(nv[1]) x + $(nv[2]) y ≥ $b")
end

function Base.in(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}
    nv, b = normal_vector(H), translation(H)
    return nv[1] * x[1] + nv[2] * x[2] ≥ b
end

function contains_in_interior(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}
    nv, b = normal_vector(H), translation(H)
    return nv[1] * x[1] + nv[2] * x[2] > b
end

Base.:(==)(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer} =
normal_vector(H1) == normal_vector(H2) && translation(H1) == translation(H2)

function Base.issubset(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer}
    normal_vector(H1) != normal_vector(H2) && return false
    return translation(H1) ≥ translation(H2)
end

@attr function direction_vector(H :: AffineHalfplane)
    nv = normal_vector(H)
    return (nv[2], -nv[1])
end

@attr function base_point(H :: AffineHalfplane)
    nv, b = normal_vector(H), translation(H)
    if nv[1] ≠ 0
        return (b // nv[1], 0 // 1)
    else
        return (0 // 1, b // nv[2])
    end
end

@attr line(H :: AffineHalfplane) = LineByDirection(base_point(H), direction_vector(H))

@attr pseudo_angle(H :: AffineHalfplane) = pseudo_angle(normal_vector(H))

@doc raw"""
    isless(H1 :: AffineHalfplane, H2 :: AffineHalfplane)

By convention, halfplanes are ordered by the angle of their normal vector
and then by distance.

"""
function Base.isless(H1 :: AffineHalfplane, H2 :: AffineHalfplane)
    isless(pseudo_angle(H1), pseudo_angle(H2)) && return true
    isless(pseudo_angle(H2), pseudo_angle(H1)) && return false
    return isless(translation(H1), translation(H2))
end



@attributes mutable struct AffineHalfplaneByNormalVector{T<:Integer} <: AffineHalfplane{T}
    normal_vector :: RationalPoint{T}
    translation :: Rational{T}

    function AffineHalfplaneByNormalVector(normal_vector :: Point{T}, translation :: Union{T, Rational{T}}) where {T <: Integer}
        !iszero(normal_vector) || error("normal vector can't be zero")
        return new{T}(normal_vector // 1, translation // 1)
    end

end

normal_vector(H :: AffineHalfplaneByNormalVector) = H.normal_vector
translation(H :: AffineHalfplaneByNormalVector) = H.translation


@attributes mutable struct AffineHalfplaneByLine{T <: Integer} <: AffineHalfplane{T}
    line :: Line{T}
    function AffineHalfplaneByLine(line :: Line{T}) where {T <: Integer}
        H = new{T}(line)
        set_attribute!(H, :line, line)
        return H
    end
end

normal_vector(H :: AffineHalfplaneByLine) = normal_vector(H.line)
direction_vector(H :: AffineHalfplaneByLine) = direction_vector(H.line)
base_point(H :: AffineHalfplaneByLine) = base_point(H.line)
line(H :: AffineHalfplaneByLine) = H.line

function translation(H :: AffineHalfplaneByLine)
    nv, p = normal_vector(H.line), base_point(H.line)
    return nv[1] * p[1] + nv[2] * p[2]
end

affine_halfplane_through_points(p :: Point{T}, q :: Point{T}) where {T <: Integer} =
AffineHalfplaneByLine(LineThroughPoints(p,q))

