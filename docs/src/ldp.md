# LDP Polygons

An LDP polygon is a lattice polygon with primitive vertices containing the
origin in its interior. LDP polygons correspond to toric log del Pezzo surfaces.
In this package, we call more generally a ``k``-rational polygon _LDP_, if it
contains the origin in its interior and its ``k``-fold multiple has primitive
vertices. With this definition, the ``k``-rational LDP polygons with exactly one
interior lattice point correspond to the toric del Pezzo surfaces with at most
``\frak{1}{k}``-log canonical singularities. Here, we list some properties of
LDP polygons that correspond to meaningful invariants of the associated toric
del Pezzo surface.

```@docs
contains_origin_in_interior
is_primitive(P :: RationalPolygon{T,N}) where {N,T <: Integer}
is_ldp
multiplicity(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
multiplicity(P :: RationalPolygon{T,N}) where {N,T <: Integer}
is_smooth
picard_index
gorenstein_index
log_canonicity
toric_prime_divisor_self_intersection
toric_prime_divisor_adjacent_intersection
degree(P :: RationalPolygon)
```

