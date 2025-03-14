# LDP polygons and toric surfaces

An LDP polygon is a lattice polygon with primitive vertices containing the
origin in its interior. LDP polygons correspond to toric _l_og _d_el _P_ezzo
surfaces. In this package, we call more generally a `k`-rational polygon LDP,
if it contains the origin in its interior and its `k`-fold dilation has
primitive vertices. Here, we list some properties of polygons that are
primarily meant to be used for LDP polygons, where they correspond to some
meaningful invariant of the associated toric del Pezzo surface.

```@docs
contains_origin_in_interior
is_primitive
is_ldp
multiplicity
is_smooth
picard_index
gorenstein_index
log_canonicity
toric_prime_divisor_self_intersection
toric_prime_divisor_adjacent_intersection
degree(P :: RationalPolygon)
```

