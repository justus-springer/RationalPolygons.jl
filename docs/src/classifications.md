# Classifications

We provide implementations for the classification algorithms of rational
polygons from Chapter ``\ref{chp:rational_polygons}``, see Sections
``\ref{doc:Maximal15090583064303023392}`` -
``\ref{doc:Almost-k-hollow-LDP-polygons}``. Additionally, we implement the
classification of LDP triangles by Picard index from Section
``\ref{sec:ldp_triangles_classification_by_picard_index}``, see
``\ref{doc:LDP-triangles-by-Picard-index}`` and the classification of LDP
quadrangles by Gorenstein index from Section
``\ref{sec:ldp_polygons_classifications_by_gorenstein_index}``, see
``\ref{doc:LDP-quadrangles-by-Gorenstein-index}``. Moreover, several other
classifications of polygons are implemented. This includes the classification of
lattice polygons by number of lattice points by Koelman [Koe91](@cite) in
Section ``\ref{doc:Lattice-polygons-by-number-of-lattice-points}``, lattice
polygons by number of interior lattice points by Castryck [Cas12](@cite) in
Section ``\ref{doc:Lattice-polygons-by-number-of-interior-lattice-points}``,
lattice polygons contained in a square by Brown and Kasprzyk [BK13](@cite) in
Section ``\ref{doc:Lattice-polygons-contained-in-a-square}``, LDP polygons by
Gorenstein index from Kasprzyk, Kreuzer and Nill [KKN10](@cite) in Section
``\ref{doc:LDP-polygons-by-Gorenstein-index}``, LDP triangles by Gorenstein
index from Bäuerle in Section ``\ref{doc:LDP-triangles-by-Gorenstein-index}``
[Ba25](@cite) and LDP triangles with integral degree by Hausen and Király
[HaKi24](@cite) in Section ``\ref{doc:LDP-triangles-with-integral-degree}``.

## Maximal rational polygons contained in ``\mathbb{R}\times[-1,1]``

We provide an implementation of Algorithm ``\ref{algo:classify_maximal_polygons_m1p1}``.

```@docs
classify_maximal_polygons_m1p1
```

## Maximal rational polygons with no interior lattice points

We provide an implementation of Algorithm ``\ref{algo:classify_maximal_polygons_no_interior_lattice_points}``.


```@docs
classify_maximal_lattice_free_polygons_m1p2_squares
classify_maximal_lattice_free_polygons_m1p2_trapezoids
classify_maximal_lattice_free_polygons_m1p2
classify_maximal_lattice_free_polygons
```

## Rational polygons with one interior lattice point

We provide an implementation of Algorithm ``\ref{algo:classify_maximal_polygons_one_interior_lattice_points}``.

```@docs
classify_maximal_polygons_genus_one_m1p1
classify_maximal_polygons_genus_one_m1p2
classify_maximal_polygons_genus_one_m2p2
classify_maximal_polygons_genus_one
classify_polygons_genus_one
```

## Almost $k$-hollow LDP polygons

We can instruct [`classify_polygons_genus_one`](@ref) to only output ``k``-rational
polygons with primitive vertices. These are exactly the almost ``k``-hollow LDP
polygons and they correspond to ``\frac{1}{k}``-log canonical toric del Pezzo surfaces.
In particular, we can reproduce the classification of the ``48032`` almost ``3``-hollow
LDP polygons (``\frac{1}{3}``-log canonical toric del Pezzo surfaces) from Theorem 4.11
of [HaHaSp25](@cite). See also Table ``\ref{class:ldp_polygons}`` for the classification up
to ``k = 6``.

```jlcon
julia> Pss = [classify_polygons_genus_one(k; primitive=true) for k = 1 : 3];

julia> numbers_of_polygons = length.(Pss)
3-element Vector{Int64}:
    16
   505
 48032

julia> max_vertices = [maximum(number_of_vertices.(Pss[k])) for k = 1 : 3]
3-element Vector{Int64}:
  6
  8
 12

julia> max_volumes = [k^2 * maximum(euclidean_area.(Pss[k])) for k = 1 : 3]
3-element Vector{Rational{Int64}}:
 9//2
 17
 47
```

## LDP triangles by Picard index

Springer [Spr24](@cite) gave an algorithm to classify LDP triangles (toric log del Pezzo
surfaces of rank one) by Picard index. `RationalPolygons.jl` implements a version
of this algorithm, which successfully reproduces the numbers from Theorem 8.5
of [Spr24](@cite).

```@docs
PicardIndexStorage
InMemoryPicardIndexStorage
HDFPicardIndexStorage
classify_lattice_triangles_by_picard_index
```

## LDP quadrangles by Gorenstein index

The following is a classification of LDP quadrangles by Gorenstein index.
A reference explaining the approach used will be added in the future.

```@docs
modified_unit_fraction_solutions
gorenstein_coefficients_to_degree_matrix_minors
classify_gorenstein_coefficients
gorenstein_coefficients_normal_form
classify_quadrilaterals_by_gorenstein_index
```

## Lattice polygons by number of lattice points

R.J. Koelman [Koe91](@cite) gave an algorithm to classify lattice polygons with
a given number of lattice points and ran it up to ``42`` lattice points, Table
4.4.3 of [Koe91](@cite). We have implemented their algorithm here, which
successfully reproduces their numbers. See also
[A371917](https://oeis.org/A371917) on OEIS for the numbers up to ``112``
lattice points.

```@docs
height_one_points
single_point_extensions
KoelmanStorage
InMemoryKoelmanStorage
HDFKoelmanStoragePreferences
HDFKoelmanStorage
classify_next_number_of_lattice_points
classify_polygons_by_number_of_lattice_points
```

## Lattice polygons by number of interior lattice points

Castryck [Cas12](@cite) gave an algorithm to classify lattice polygons by number
of interior lattice points and ran it up to 30 interior lattice points. Our
implementation here successfully reproduces their numbers from Table 1 of
[Cas12](@cite), see also [A322343](https://oeis.org/A322343) on OEIS.

```@docs
classify_maximal_lattice_polygons_with_collinear_interior_points
classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior
CastryckStorage
InMemoryCastryckStorage
HDFCastryckStoragePreferences
HDFCastryckStorage
classify_next_genus
classify_lattice_polygons_by_genus
```

## Lattice polygons contained in a square

Brown and Kasprzyk [BK13](@cite) authors considered lattice polygons that are
contained in a square of fixed side length and classified them up to side length
``7``. Their numbers (Table 1 of [BK13](@cite), see also
[A374975](https://oeis.org/A374975)) can be reproduced with
`RationalPolygons.jl` as follows:

```julia
julia> square(m) = convex_hull(LatticePoint{Int}[(0,0),(m,0),(0,m),(m,m)])
square (generic function with 1 method)

julia> Pss = [subpolygons(square(m); use_affine_normal_form = true, only_equal_number_of_interior_lattice_points = false) for m = 1 : 7];

julia> numbers_of_polygons = [length(Pss[m]) - length(Pss[m-1]) for m = 2 : 7]
6-element Vector{Int64}:
      15
     131
    1369
   13842
  129185
 1104895

julia> max_vertices = [maximum(number_of_vertices.(Pss[m])) for m = 1 : 7]
7-element Vector{Int64}:
  4
  6
  8
  9
 10
 12
 13

julia> [length(filter(P -> number_of_vertices(P) == max_vertices[m], Pss[m])) for m = 1 : 7]
7-element Vector{Int64}:
  1
  1
  1
  1
 15
  2
  3
```

## LDP polygons by Gorenstein index

Kasprzyk, Kreuzer and Nill [KKN10](@cite) gave an algorithm to classify LDP
polygons by Gorenstein index. `RationalPolygons.jl` implements a version of their
algorithm, which successfully reproduces their numbers (see Theorem 1.2 of
[KKN10](@cite)).

```@docs
PartialLDP
initial_special_facets
choose_next_vertex
classify_lattice_polygons_by_gorenstein_index
```

## LDP triangles by Gorenstein index

Bäuerle [Ba25](@cite) gave an algorithm to classify Fano simplices by dimension
and Gorenstein index. `RationalPolygons.jl` implements a version of his
algorithm (specialized to the two-dimensional case), which reproduces his
numbers successfully (see Theorem 1.4 of [Ba25](@cite) and
[A145582](https://oeis.org/A145582)).

```@docs
unit_fraction_partitions_length_three
BaeuerleStorage
InMemoryBaeuerleStorage
HDFBaeuerleStorage
classify_lattice_triangles_by_gorenstein_index
```

## LDP triangles with integral degree

Hausen and Király [HaKi24](@cite) classified fake weighted projective planes
having integral degree (=canonical self intersection). In terms of polygons,
these can be described as LDP triangles such that twice the euclidean area of
its dual is an integer. The attained values of this integer (which is the degree
of the associated fake weighted projective plane) are the integers from ``1`` to
``9``, excluding ``7``. In total, there are 24 infinite series of these
triangles, where each of them is parameterized by the solution set of a squared
Markov type equation (see Theorem 1.1 of [HaKi24](@cite)). These solution sets
can be described as infinite binary trees with a unique root.
`RationalPolygons.jl` uses this description to implement a classification
algorithm for LDP triangles with integral degree. To make this classification
finite, one has to provide a maximal depth to which the solution trees are
traversed.

All methods used for this classification come with a parameter `T <: Integer`,
which is the integer type to be used. Since the entries of the solution triples
of squared Markov type equations grow very quickly with increasing depth, it is
recommended to use `BigInt` here instead of fixed-size integer types (`Int64`
overflows already for `depth = 4`).

```@docs
degree(w :: SVector{3})
fake_weight_vector
n_step_mutations
initial_triple
adjust_triple
classify_squared_markov_type_equation_solutions
fake_weight_vectors_to_triangles
classify_lattice_triangles_integral_degree
```

