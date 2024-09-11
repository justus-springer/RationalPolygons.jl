# 2D Geometry

RationalPolygons.jl comes with its own library for two-dimensional geometry over
the rational numbers, which is implented from scratch in pure Julia.

## Points

The basic types for points in RationalPolygons.jl are `LatticePoint`, `RationalPoint` and `Point`, which are aliases for [static vectors](https://juliaarrays.github.io/StaticArrays.jl/stable/) of length two.

```@docs
LatticePoint
RationalPoint
Point
is_k_rational(k :: T, p :: Point{T}) where {T <: Integer}
is_integral
rationality(p :: Point)
multiplicity
is_primitive(p :: Point)
norm
distance
pseudo_angle
```

## Graham scan


```@docs
graham_scan!
graham_scan
```

## Lines

```@docs
Line
base_point
direction_vector
line_through_points
horizontal_line
vertical_line
Base.in(x :: Point{T}, L :: Line{T}) where {T <: Integer}
normal_vector(L :: Line{T}) where {T <: Integer}
IntersectionBehaviour
IntersectInPoint
NoIntersection
LinesAreEqual
intersection_behaviour
intersection_point
```

## Affine halfplanes

```@docs
AffineHalfplane
affine_halfplane
normal_vector(H :: AffineHalfplane{T}) where {T <: Integer}
translation
Base.in(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}
contains_in_interior(x :: Point{T}, H :: AffineHalfplane{T}) where {T <: Integer}
Base.issubset(H1 :: AffineHalfplane{T}, H2 :: AffineHalfplane{T}) where {T <: Integer}
line(H :: AffineHalfplane{T}) where {T <: Integer}
```


