
@doc raw"""
    abstract type AffineHalfplane end

An affine halfplane in two-dimensional rational space. Subtypes must at least
implement `normal_vector` and `translation`.

"""
struct AffineHalfplane{T <: Integer}
    normal_vector :: RationalPoint{T}
    translation :: Rational{T}

    function AffineHalfplane(normal_vector :: Point{T}, translation :: Union{T, Rational{T}}) where {T <: Integer}
        !iszero(normal_vector) || error("normal vector can't be zero")
        return new{T}(normal_vector, translation)
    end

end

affine_halfplane(normal_vector :: Point{T}, translation :: Union{T, Rational{T}}) where {T <: Integer} =
AffineHalfplane(normal_vector, translation)

function affine_halfplane(L :: Line{T}) where {T <: Integer}
    nv, p = normal_vector(L), base_point(L)
    return affine_halfplane(nv, nv[1] * p[1] + nv[2] * p[2])
end

affine_halfplane(p :: Point{T}, q :: Point{T}) where {T <: Integer} =
affine_halfplane(line_through_points(p,q))

normal_vector(H :: AffineHalfplane{T}) where {T <: Integer} = H.normal_vector
translation(H :: AffineHalfplane{T}) where {T <: Integer} = H.translation

function Base.hash(H :: AffineHalfplane{T}, h :: UInt64) where {T <: Integer}
    h = hash(normal_vector(H), h)
    h = hash(translation(H), h)
    return h
end

Base.:(==)(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer} =
normal_vector(H1) == normal_vector(H2) && translation(H1) == translation(H2)

function Base.show(io :: IO, H :: AffineHalfplane{T}) where {T <: Integer}
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

function Base.issubset(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer}
    normal_vector(H1) != normal_vector(H2) && return false
    return translation(H1) ≥ translation(H2)
end

function direction_vector(H :: AffineHalfplane{T}) where {T <: Integer}
    nv = normal_vector(H)
    return RationalPoint(nv[2], -nv[1])
end

function base_point(H :: AffineHalfplane{T}) where {T <: Integer}
    nv, b = normal_vector(H), translation(H)
    if nv[1] ≠ 0
        return RationalPoint{T}(b // nv[1], T(0) // T(1))
    else
        return RationalPoint{T}(T(0) // T(1), b // nv[2])
    end
end

line(H :: AffineHalfplane{T}) where {T <: Integer} = Line(base_point(H), direction_vector(H))

pseudo_angle(H :: AffineHalfplane{T}) where {T <: Integer} = pseudo_angle(normal_vector(H))

@doc raw"""
    isless(H1 :: AffineHalfplane, H2 :: AffineHalfplane)

By convention, halfplanes are ordered by the angle of their normal vector
and then by distance.

"""
function Base.isless(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer}
    isless(pseudo_angle(H1), pseudo_angle(H2)) && return true
    isless(pseudo_angle(H2), pseudo_angle(H1)) && return false
    return isless(translation(H1), translation(H2))
end

function rand(::Type{AffineHalfplane}, nv_range, translation_range)
    nv = (rand(nv_range),rand(nv_range))
    while iszero(nv)
        nv = (rand(nv_range),rand(nv_range))
    end
    nv = primitivize(nv)
    t = rand(translation_range)
    return AffineHalfplaneByNormalVector(nv, t)
end

Base.:(+)(H :: AffineHalfplane{T}, x :: Union{T, Rational{T}}) where {T <: Integer} =
affine_halfplane(normal_vector(H), translation(H) + x)

Base.:(-)(H :: AffineHalfplane{T}, x :: Union{T, Rational{T}}) where {T <: Integer} =
affine_halfplane(normal_vector(H), translation(H) - x)

Base.:(-)(H :: AffineHalfplane{T}) where {T <: Integer} =
affine_halfplane(-normal_vector(H), -translation(H))
