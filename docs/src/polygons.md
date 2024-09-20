# Rational Polygons

## The `RationalPolygon` type

In RationalPolygons.jl, we represent a polygon ``P \subseteq \mathbb{R}^2`` by
two pieces of data: An integral matrix ``V \in \mathbb{Z}^{2\times N}`` and an
integer ``k \in \mathbb{Z}``, called the _rationality_. The associated polygon
has as vertices the columns of ``V`` divided by ``k``. To represent ``V``, we
use a [static matrix](https://juliaarrays.github.io/StaticArrays.jl/stable/),
which are faster than Julia's internal matrices for many common operations. 

The standard lattice triangle can be created as follows:

```julia
julia> using RationalPolygons, StaticArrays

julia> P = RationalPolygon(SMatrix{2,3}([0 1 0 ; 0 0 1]), 1)
Rational polygon of rationality 1 with 3 vertices.
```

!!! warning 
    When creating a `RationalPolygon` from a constructor, the user has
    to be certain that the columns of ``V`` truly are vertices of a polygons, i.e.
    no column is contained in the convex hull of the other columns and the columns
    are sorted by angle (both clockwise and counterclockwise is allowed). If this
    is not known, use [`convex_hull`](@ref) to create the polygon instead.

There are two ways in which this encoding is not unique: First, scaling ``V``
and ``k`` by the same factor does not change the polygon, e.g. ``(V,k)`` describes
the same polygon as ``(2V,2k)``. Even though they are
mathematically the same polygon, RationalPolygon.jl views them as different
objects, once viewed as a ``k``-rational polygon and once viewed as a
``2k``-rational polygon. The second way in which this encoding is not unique is
that there is no canonical "first vertex" of a polygon, i.e. we can shift the
columns of ``V`` around and still describe the same polygon. Moreover, we
choose to order them clockwise or counterclockwise. This problem is adressed in
the section on [normal forms](##Normal forms).

### Constructors

Besides the type constructor methods, we provide the functions
[`convex_hull`](@ref) and [`intersect_halfplanes`](@ref) to create a polygon
from an unstructured collection of points in the plane of affine halfplanes.

```@docs
RationalPolygon
convex_hull
intersect_halfplanes
empty_polygon
```

### Properties

```@docs
number_of_vertices
rationality(P :: RationalPolygon)
vertex_matrix
is_unimodular_normal_form
is_affine_normal_form
scaled_vertex
vertex
vertices
affine_halfplane(P :: RationalPolygon, i :: Int)
affine_halfplanes(P :: RationalPolygon)
Base.in(x :: Point{T}, P :: RationalPolygon{T}) where {T <: Integer}
contains_in_interior(x :: Point{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer}
dim
is_primitive(P :: RationalPolygon{T,N}) where {N,T <: Integer}
is_fano(P :: RationalPolygon)
normalized_area
euclidian_area
is_maximal
dual
gorenstein_index
log_canonicity
```

## Ehrhart Theory

Consider a ``k``-rational polygon ``P``. The main result of Ehrhart Theory is
that the the number of lattice points in integral multiples of ``P`` is a
quasipolynomial, called its _Ehrhart quasipolynomial_:

```math
\mathrm{ehr_P}(t) = |tP \cap \mathbb{Z}^2| = At^2 + a(t)t+b(t), \qquad t \in \mathbb{Z}.
```

Here, ``A`` is the euclidian area of ``P`` and ``a, b\colon \mathbb{Z} \to
\mathbb{Q}`` are ``k``-periodic functions. These can be computed by

```math
a(t) = -(2t+k)\cdot A + \frac{\mathrm{ehr}_P(t+k)-\mathrm{ehr}_P(t)}{k},
```
```math
b(t) = (t^2+tk)\cdot A + \frac{(t+k)\mathrm{ehr}_P(t)-t\mathrm{ehr}_P(t+k)}{k}.
```

Setting ``\tilde{A} := 2k^2A,\ \tilde{a} := 2k^2 a`` and ``\tilde{b} := 2k^2b``, we get
integer valued functions ``\tilde a`` and ``\tilde b`` which we call the
_normalized Ehrhart coefficients_. We then encode the Ehrhart quasipolynomial
by the ``3\times k`` integral matrix of its normalized Ehrhart coefficients:

```math
\begin{bmatrix}
\tilde{A} & \tilde{a}(1) & \tilde{b}(1) \\
\tilde{A} & \tilde{a}(2) & \tilde{b}(2) \\
\vdots & \vdots & \vdots \\
\tilde{A} & \tilde{a}(k) & \tilde{b}(k)
\end{bmatrix} \in \mathbb{Z}^{3\times k}.
```

If ``P`` is integral, we have ``k=1`` and its Ehrhart quasipolynomial is a
regular polyomial of degree 2. In general, the periods of ``a`` and ``b`` are
divisors of ``k``. If they are strictly smaller than ``k``, we speak of _quasiperiod collapse_. A rational polygon is called _quasiintegral_ if the
periods of ``a`` and ``b`` are both 1, hence it has an Ehrhart polynomial.

RationalPolygons.jl comes with many methods for counting the (interior,
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


## Normal forms and automorphism groups

Two ``k``-rational polygons are called _(affine) unimodular_ equivalent if they
can be transformed into each other by an (affine) unimodular transformation.
The purpose of a normal form is to provide a unique representative for every
equivalence class, i.e. two polygons should be (affine) unimodular equivalent
to each other if and only if their (affine) unimodular normal forms coincide.

For details about the normal form used in RationalPolygons.jl, we refer to
[BS24](@cite).

### Normal forms

```@docs
unimodular_normal_form
are_unimodular_equivalent
affine_normal_form
are_affine_equivalent
```

### Automorphism groups

```@docs
PolygonAutomorphismGroup
CyclicGroup
DihedralGroup
is_cyclic
order
unimodular_automorphism_group
affine_automorphism_group
```

## Lattice width

In [Boh23](@cite), Bohnert describes the concept of _lattice width data_, which
captures information about the slicing lengths of a polygon with respect to a
given direction vectors. RationalPolygons.jl implements this concept, following his Definition 2.11.

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

RationalPolygons.jl provides two ways to save polygons to a file: The first is text-based, where polygons can be written and read to files containing one polygon per line like this:

```shell
[[2, 0], [1, 3], [-1, 0], [-3, -4]]
[[1, 0], [2, 6], [-4, -9]]
[[1, 0], [3, 5], [0, 1], [-5, -8]]
....
```

This text-based format has the advantage of being universally understandable
and easy to use. However, storing polygons as ascii strings is not very
space-efficients, as they contain lots of redundant control characters. Hence
we provide another way to store polygons in binary form, which uses the HDF5
format and is more suitable for large datasets. For an example session, see
[`write_polygon_dataset`](@ref).

```@docs
parse_rational_polygons
write_rational_polygons
create_polygon_dataset
write_polygon_dataset
read_polygon_dataset
```

