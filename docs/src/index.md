This appendix serves as a reference for `RationalPolygons.jl`, a Julia package
for computations with rational convex polygons. `RationalPolygons.jl` does not
make use of any external computer algebra system but implements all necessary
algorithms, including two-dimensional euclidean geometry, from scratch in pure
Julia. This allows for quite good performance, with computations involving
billions of polygons being feasible on a personal computer.

Let us give an impression of `RationalPolygons.jl` by means of an example session.
After loading in the package, we can create a polygon by calling the
[`convex_hull`](@ref) function on a set of rational points:

```@repl quick_start
using RationalPolygons
P = convex_hull(RationalPoint{Int}[(-3//2,-1//2), (-1,-1), (1//2,-1//2), (3//2,1//2), (1,1), (-1//2,1//2)])
```

Using Julia's plotting library, we can visualize the polygon:

```@repl quick_start
using Plots
plot(P);
```

![image](example_polygon.png)

Next, we compute some basic properties:

```@repl quick_start
number_of_interior_lattice_points(P)
number_of_boundary_lattice_points(P)
euclidean_area(P)
ehrhart_quasipolynomial(P)
affine_automorphism_group(P)
gorenstein_index(P)
```

Lastly, we verify that the polygon is equivalent to its own dual.

```@repl quick_start
dual(P)
are_affine_equivalent(P, dual(P))
```

The rest of this appendix is organized as follows: In Section
``\ref{doc:2D-Geometry}``, we go over some basic functions for two-dimensional
euclidean geometry over the rationals, such as computing the convex hull and
intersecting lines. Section ``\ref{doc:Polygons}`` covers the type of rational
polygons as well as basic properties and the normal form. Section
``\ref{doc:LDP-Polygons}`` is about LDP polygons and their relation to toric del
Pezzo surfaces. In Section ``\ref{doc:Subpolygons}``, we discuss computation
subpolygons, following the approach from Section ``\ref{subsec:subpolygons}``.
Finally, Section ``\ref{doc:Classifications}`` covers implementations of the
classification algorithms from Chapter ``\ref{chp:rational_polygons}``, Sections
``\ref{sec:ldp_triangles_classification_by_picard_index}`` and
``\ref{sec:ldp_polygons_classifications_by_gorenstein_index}``, as well as
various classification algorithms by other authors.

This documentation has been generated from the docstrings of the package's
source code. A web version is available on its GitHub page
[RationalPolygons_jl](@cite). All example sessions have been tested against
version `v1.2.0`.
