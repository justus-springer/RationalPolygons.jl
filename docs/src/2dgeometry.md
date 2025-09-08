# 2D Geometry

We provide basic functionality for two-dimensional geometry over the rational
numbers. This includes distances and angles, Graham scan for computing the
convex hull, intersection of lines, and affine halfplanes.

## Points

In `RationalPolygons.jl`, we represent points as [static
vectors](https://juliaarrays.github.io/StaticArrays.jl/stable/) of length two.
As the name suggests, these are _statically sized_, which leads to improved
performance and memory management for many common operations. We provide three
type aliases [`LatticePoint`](@ref), [`RationalPoint`](@ref), and
[`Point`](@ref), the latter being the union of the previous two. Note that
everything is stated for an arbitrary subtype `T <: Integer`, e.g. fixed size
machine integers like `Int64`, or unbounded integer types like `BigInt`.

```@docs
LatticePoint
RationalPoint
Point
is_k_rational(k :: T, p :: Point{T}) where {T <: Integer}
is_integral
Base.denominator(p :: Point)
is_primitive(p :: Point)
multiplicity(p :: Point)
primitivize(p :: Point)
norm
distance
pseudo_angle(p :: Point{T}) where {T <: Integer}
```

## Graham scan

The Graham scan is a planar convex hull algorithm named after Ronald Graham
[Gra72](@cite). With an asymptotic running time of ``O(n \cdot \mathrm{log}(n))``, it
is a lot quicker than algorithms that work in arbitrary dimension.

```@docs
graham_scan!
graham_scan
```

## Lines

We provide basic functionality for lines in the two-dimensional rational plane.
A line is encoded as a simple struct consisting of a base point and a direction
vector. Two lines can be intersected to produce a value of type
[`IntersectionBehaviour`](@ref). This is an abstract type with three concrete
subtypes [`NoIntersection`](@ref), [`LinesAreEqual`](@ref) and
[`IntersectInPoint`](@ref), capturing the three intersection scenarios.

```@docs
Line
base_point(L :: Line{T}) where {T <: Integer}
direction_vector(L :: Line{T}) where {T <: Integer}
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

An affine halfplane is encoded as a struct consisting of a normal vector ``v \in
\mathbb{Q}^2`` and a translation ``b \in \mathbb{Q}^2``, representing the set of
points ``\{ x \in \mathbb{Q}^2 | \langle v, x \rangle \geq b \}``.

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


