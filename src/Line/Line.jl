
@doc raw"""
    Line{T <: Integer}

A line in 2-dimensional rational space. It has two fields: `base_point :: RationalPoint{T}` and `direction_vector :: RationalPoint{T}`.

"""
struct Line{T <: Integer}
    base_point :: RationalPoint{T}
    direction_vector :: RationalPoint{T}
    function Line(base_point :: Point{T}, direction_vector :: Point{T}) where {T <: Integer}
        !iszero(direction_vector) || error("direction vector can't be zero")
        new{T}(base_point, direction_vector)
    end
end


@doc raw"""
    base_point(L :: Line{T}) where {T <: Integer}

Return a point on the line `L`.

"""
base_point(L :: Line{T}) where {T <: Integer} = L.base_point


@doc raw"""
    direction_vector(L :: Line{T}) where {T <: Integer}

Return the direction vector of `L`.

"""
direction_vector(L :: Line{T}) where {T <: Integer} = L.direction_vector

Base.show(io :: IO, L :: Line) =
print(io, "Line with base point $(base_point(L)) and direction vector $(direction_vector(L))")

Base.:(==)(L1 :: Line, L2 :: Line) =
base_point(L1) == base_point(L2) && direction_vector(L1) == direction_vector(L2)

function Base.hash(L :: Line, h :: UInt64) 
    h = hash(L.base_point, h)
    h = hash(L.direction_vector, h)
    return h
end


@doc raw"""
    Base.in(x :: Point{T}, L :: Line{T}) where {T <: Integer}

Check whether a point `x` lies on a line `L`.

# Example

```jldoctest
julia> RationalPoint(3//2,1) ∈ line_through_points(Point(1,0),Point(2,2))
true
```

"""
Base.in(x :: Point{T}, L :: Line{T}) where {T <: Integer} =
iszero(det(x - base_point(L), direction_vector(L)))


@doc raw"""
    normal_vector(L :: Line{T}) where {T <: Integer}

Return a primitive vector orthogonal to the direction vector of `L`.

# Example

```jldoctest
julia> normal_vector(line_through_points(Point(1,0),Point(2,2)))
2-element StaticArraysCore.SVector{2, Int64} with indices SOneTo(2):
 -2
  1
```

"""
function normal_vector(L :: Line{T}) where {T <: Integer}
    v = primitivize(direction_vector(L))
    return LatticePoint(-v[2], v[1])
end


@doc raw"""
    line_through_points(A :: Point{T}, B :: Point{T}) where {T <: Integer}

Return the line going through the points `A` and `B`.

"""
line_through_points(A :: Point{T}, B :: Point{T}) where {T <: Integer} =
Line(A, B - A)


@doc raw"""
    horizontal_line(y :: Union{T, Rational{T}}) where {T <: Integer}

Return the horizontal line at `y`.

"""
horizontal_line(y :: Union{T, Rational{T}}) where {T <: Integer} =
Line(Point(zero(y), y), Point(one(y), zero(y)))


@doc raw"""
    vertical_line(x :: Union{T, Rational{T}}) where {T <: Integer}

Return the vertical line at `x`.

"""
vertical_line(x :: Union{T, Rational{T}}) where {T <: Integer} =
Line(Point(x, zero(x)), Point(zero(x), one(x)))
