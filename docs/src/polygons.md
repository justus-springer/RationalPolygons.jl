# Polygons

In `RationalPolygons.jl`, we represent a polygon ``P \subseteq \mathbb{R}^2`` by
two pieces of data: An integral matrix ``V \in \mathbb{Z}^{2\times N}``, called
the _vertex matrix_ and an integer ``k \in \mathbb{Z}``, called the
_rationality_. The associated polygon has as vertices the columns of ``V``
divided by ``k``. To represent ``V``, we use a [static
matrix](https://juliaarrays.github.io/StaticArrays.jl/stable/), which are faster
than Julia's internal matrices for many common operations. However, this implies
that the type `RationalPolygon` will depend on the number of vertices `N` as a
type parameter, which affects Julia's dispatch mechanism at runtime. As long as
the number of distinct values of `N` occurring during a computation remains
relatively small, this should not cause performance issues.

There are two ways in which our encoding of rational polygons is not unique:
First, scaling ``V`` and ``k`` by the same factor does not change the polygon,
e.g. ``(V,k)`` describes the same polygon as ``(2V,2k)``. Even though they are
mathematically the same polygon, `RationalPolygon.jl` views them as different
objects: once as a ``k``-rational polygon and once as a ``2k``-rational polygon.
The second way in which this encoding is not unique is that we can change the
order of the columns. While we require them to be sorted counterclockwise, we
may use any vertex as the first column.

```@docs
RationalPolygon
```

The rest of this section is organized as follows. In Subsection
``\ref{doc:Constructors}``, we describe different methods of constructing a
rational polygon. Subsection ``\ref{doc:Basic-Properties}`` is about basic
properties. In Subsection ``\ref{doc:Ehrhart-Theory}``, we discuss various ways
of counting lattice points as well as Ehrhart theory. Subsection
``\ref{doc:LDP-polygons-and-toric-surfaces}`` is about properties of LDP
polygons that correspond to meaningful invariants of their associated toric del
Pezzo surfaces. In Subsection
``\ref{doc:Normal-forms-and-automorphism-groups}``, we discuss unimodular and
affine unimodular normal forms. Finally, Subsection ``\ref{doc:Lattice-width}``
is about lattice width and related concepts.

## Constructors

To construct a rational polygon, one can either use a type constructor
or one of the functions [`convex_hull`](@ref) and
[`intersect_halfplanes`](@ref). We describe the type constructors first.

```@docs
RationalPolygon(vertex_matrix :: SMatrix{2,N,T,M}, rationality :: T) where {N, M, T <: Integer}
convex_hull
intersect_halfplanes
empty_polygon
```

## Basic Properties

We provide basic properties and checks for rational polygons. Note that indices
corresponding to vertices are always considered cyclic, i.e. the ``N+1``-th vertex
of a polygon with ``N`` vertices cycles back to its first vertex.

```@docs
number_of_vertices
rationality(P :: RationalPolygon)
Base.denominator(P :: RationalPolygon)
vertex_matrix
scaled_vertex
vertex
vertices
affine_halfplane(P :: RationalPolygon, i :: Int)
affine_halfplanes(P :: RationalPolygon)
Base.in(x :: Point{T}, P :: RationalPolygon{T}) where {T <: Integer}
contains_in_interior(x :: Point{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer}
dim
normalized_area
euclidean_area
is_maximal
dual
```


## Ehrhart Theory

Recall that the number of lattice points in integral multiples of a
``k``-rational polygon ``P`` is a quasipolynomial, called its _Ehrhart
quasipolynomial_:

```math
\mathrm{ehr_P}(t) = |tP \cap \mathbb{Z}^2| = At^2 + a(t)t+b(t), \qquad t \in \mathbb{Z}.
```

Here, ``A`` is the euclidean area of ``P`` and ``a, b\colon \mathbb{Z} \to
\mathbb{Q}`` are ``k``-periodic functions. These can be computed by

```math
a(t) = -(2t+k)\cdot A + \frac{\mathrm{ehr}_P(t+k)-\mathrm{ehr}_P(t)}{k},
```
```math
b(t) = (t^2+tk)\cdot A + \frac{(t+k)\mathrm{ehr}_P(t)-t\mathrm{ehr}_P(t+k)}{k}.
```

Setting ``\tilde{A} := 2k^2A,\ \tilde{a} := 2k^2 a`` and ``\tilde{b} := 2k^2b``,
we get integer-valued functions ``\tilde a`` and ``\tilde b``, which we call the
_normalized Ehrhart coefficients_. In `RationalPolygons.jl`, we encode the
Ehrhart quasipolynomial by the ``3\times k``-integral matrix of its normalized
Ehrhart coefficients:

```math
\begin{bmatrix}
\tilde{A} & \tilde{a}(1) & \tilde{b}(1) \\
\tilde{A} & \tilde{a}(2) & \tilde{b}(2) \\
\vdots & \vdots & \vdots \\
\tilde{A} & \tilde{a}(k) & \tilde{b}(k)
\end{bmatrix}
```

If ``P`` is integral, we have ``k=1`` and its Ehrhart quasipolynomial is a
regular polynomial of degree 2. In general, the periods of ``a`` and ``b`` are
divisors of ``k``. If they are strictly smaller than ``k``, we speak of
_quasiperiod collapse_. A rational polygon is called _quasiintegral_ if the
periods of ``a`` and ``b`` are both 1, hence it has an Ehrhart polynomial.

`RationalPolygons.jl` comes with many methods for counting the (interior,
boundary) lattice points of a rational polygon as well as computing its
Ehrhart quasipolynomial and its periods.

### Counting lattice points

```@docs
generic_lattice_points
boundary_k_rational_points
number_of_boundary_k_rational_points
boundary_lattice_points
number_of_boundary_lattice_points
interior_k_rational_points
number_of_interior_k_rational_points
interior_lattice_points
number_of_interior_lattice_points
k_rational_points
number_of_k_rational_points
lattice_points
number_of_lattice_points
k_rational_hull
interior_k_rational_hull
integer_hull
interior_integer_hull
```

### Ehrhart quasipolynomial

```@docs
is_periodic
period
ehrhart_quasipolynomial_with_periods
ehrhart_quasipolynomial
ehrhart_quasipolynomial_periods
ehrhart_quasipolynomial_period
is_quasiintegral
```

## LDP polygons and toric surfaces

An LDP polygon is a lattice polygon with primitive vertices containing the
origin in its interior. LDP polygons correspond to toric log del Pezzo surfaces.
In this package, we call more generally a ``k``-rational polygon _LDP_, if it
contains the origin in its interior and its ``k``-fold multiple has primitive
vertices. With this definition, the ``k``-rational LDP polygons with exactly one
interior lattice point correspond to the toric del Pezzo surfaces with at most
``\frac{1}{k}``-log canonical singularities. Here, we list some properties of
LDP polygons that correspond to meaningful invariants of the associated toric
del Pezzo surface. For more background on LDP polygons, see Section ``\ref{sec:ldp_polygons_toric_log_del_pezzo_surfaces}``.

```@docs
contains_origin_in_interior
is_primitive(P :: RationalPolygon{T,N}) where {N,T <: Integer}
is_ldp
multiplicity(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
multiplicity(P :: RationalPolygon{T,N}) where {N,T <: Integer}
grading_matrix_free_part
grading_matrix_torsion_part
grading_matrix
is_smooth
picard_index
gorenstein_index
gorenstein_coefficients
gorenstein_matrix
log_canonicity
toric_prime_divisor_self_intersection
toric_prime_divisor_adjacent_intersection
degree(P :: RationalPolygon)
```


## Normal forms and automorphism groups

Two ``k``-rational polygons are called _(affine) unimodular_ equivalent if they
can be transformed into each other by an (affine) unimodular transformation.
The purpose of a normal form is to provide a unique representative for every
equivalence class, i.e. two polygons should be (affine) unimodular equivalent
to each other if and only if their (affine) unimodular normal forms coincide.

For details about the normal form used in `RationalPolygons.jl`, see Section
``\ref{subsec:normal_forms}``.

```@docs
unimodular_normal_form
are_unimodular_equivalent
affine_normal_form
are_affine_equivalent
PolygonAutomorphismGroup
CyclicGroup
DihedralGroup
is_cyclic
order
unimodular_automorphism_group
affine_automorphism_group
```

## Lattice width

We provide functions to compute the lattice width as well as all direction vectors
in which the lattice width is attained. Furthermore, we implement the concept of
_lattice width data_ following [Boh23](@cite), which captures information about
the slicing lengths of a polygon with respect to given direction vectors.

```@docs
width
all_direction_vectors_with_width_less_than
width_direction_vectors
adjust_to_width_direction
number_of_interior_integral_lines
minimal_number_of_interior_integral_lines
is_realizable_in_interval
LatticeWidthData
lattice_width_data
number_of_interior_integral_vertical_lines
position_of_longest_vertical_slice_length
lattice_width_datas
numbers_of_interior_integral_vertical_lines
positions_of_longest_vertical_slice_length
```

## IO

`RationalPolygons.jl` provides two ways to read and write polygons from files.
The first is text-based. Polygons can be read from files containing one polygon
per line like this:

```shell
[[2, 0], [1, 3], [-1, 0], [-3, -4]]
[[1, 0], [2, 6], [-4, -9]]
[[1, 0], [3, 5], [0, 1], [-5, -8]]
....
```

This text-based format has the advantage of being easy to understand and use.
However, storing polygons as ASCII strings is not very space-efficient, as they
contain lots of redundant control characters. Hence we provide another way to
store polygons in binary and compressed form, which uses the HDF5 file format.
This is more suitable for large datasets. For an example session, see
[`write_polygon_dataset`](@ref).

```@docs
parse_rational_polygons
write_rational_polygons
create_polygon_dataset
write_polygon_dataset
read_polygon_dataset
```

