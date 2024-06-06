
@doc raw"""
    AffineHalfplane{T <: Integer}

An affine halfplane in two-dimensional rational space.

"""
struct AffineHalfplane{T <: Integer}
    normal_vector :: RationalPoint{T}
    translation :: Rational{T}

    function AffineHalfplane(normal_vector :: Point{T}, translation :: Union{T, Rational{T}}) where {T <: Integer}
        !iszero(normal_vector) || error("normal vector can't be zero")
        return new{T}(normal_vector, translation)
    end

end


@doc raw"""
    affine_halfplane(nv :: Point{T}, b :: Union{T, Rational{T}}) where {T <: Integer}

Return the affine halfplane given by the equation `nv[1] * x[1] + nv[2] * x[2]
≥ b`.

"""
affine_halfplane(nv :: Point{T}, b :: Union{T, Rational{T}}) where {T <: Integer} =
AffineHalfplane(nv, b)


@doc raw"""
    affine_halfplane(L :: Line{T}) where {T <: Integer}

Return the affine halfplane associated to a line in 2D space. The halfplane is
understood to consist of those points *to the left* of the line `L`, looking in
the direction given by `direction_vector(L)`.

"""
function affine_halfplane(L :: Line{T}) where {T <: Integer}
    nv, p = normal_vector(L), base_point(L)
    return affine_halfplane(nv, nv[1] * p[1] + nv[2] * p[2])
end


@doc raw"""
    affine_halfplane(p :: Point{T}, q :: Point{T}) where {T <: Integer}

Return the affine halfplane associated to the line going through the points `p`
and `q`.

"""
affine_halfplane(p :: Point{T}, q :: Point{T}) where {T <: Integer} =
affine_halfplane(line_through_points(p,q))


@doc raw"""
    normal_vector(H :: AffineHalfplane{T}) where {T <: Integer}

Return the normal vector of the `H`.

"""
normal_vector(H :: AffineHalfplane{T}) where {T <: Integer} = H.normal_vector


@doc raw"""
    translation(H :: AffineHalfplane{T}) where {T <: Integer}

Return the affine translation of `H`.

"""
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


@doc raw"""
    Base.in(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}

Check whether a point `x` lies in the halfplane `H`.

"""
function Base.in(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}
    nv, b = normal_vector(H), translation(H)
    return nv[1] * x[1] + nv[2] * x[2] ≥ b
end


@doc raw"""
    contains_in_interior(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}

Check whether a point `x` lies in the interior of `H`.

"""
function contains_in_interior(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}
    nv, b = normal_vector(H), translation(H)
    return nv[1] * x[1] + nv[2] * x[2] > b
end


@doc raw"""
    Base.issubset(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer}

Check whether `H1` is a subset of `H2`.

"""
function Base.issubset(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer}
    normal_vector(H1) != normal_vector(H2) && return false
    return translation(H1) ≥ translation(H2)
end


@doc raw"""
    direction_vector(H :: AffineHalfplane{T}) where {T <: Integer}

Return the direction vector of the line associated to `H`.

"""
function direction_vector(H :: AffineHalfplane{T}) where {T <: Integer}
    nv = normal_vector(H)
    return RationalPoint(nv[2], -nv[1])
end


@doc raw"""
    base_point(H :: AffineHalfplane{T}) where {T <: Integer}

Return a point on the line associated to `H`.

"""
function base_point(H :: AffineHalfplane{T}) where {T <: Integer}
    nv, b = normal_vector(H), translation(H)
    if nv[1] ≠ 0
        return RationalPoint{T}(b // nv[1], T(0) // T(1))
    else
        return RationalPoint{T}(T(0) // T(1), b // nv[2])
    end
end


@doc raw"""
    line(H :: AffineHalfplane{T}) where {T <: Integer}

Return the line associated to `H`.

"""
line(H :: AffineHalfplane{T}) where {T <: Integer} = Line(base_point(H), direction_vector(H))


@doc raw"""
    pseudo_angle(H :: AffineHalfplane{T}) where {T <: Integer}

Return the pseudo angle of the normal vector of `H`.

"""
pseudo_angle(H :: AffineHalfplane{T}) where {T <: Integer} = pseudo_angle(normal_vector(H))


@doc raw"""
    isless(H1 :: AffineHalfplane, H2 :: AffineHalfplane)

Test whether `H1` comes before `H2`, where we order them first by the angle of
their normal vector and then by their affine translation.

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
