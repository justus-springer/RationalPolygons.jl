# Subpolygons

A polygon ``P`` is called a _subpolygon_ of ``Q``, if ``\varphi(P) \subseteq Q``
for some affine unimodular transformation ``\varphi``. Given a ``k``-rational
polygon ``P``, we can find all ``k``-rational subpolygons of ``P`` using the algorithm from
Subsection ``\ref{subsec:subpolygons}``. The idea is to successively remove
vertices of ``P`` by computing Hilbert bases.

## Hilbert bases

We follow [CoLiSc11](@cite) to compute Hilbert bases of two-dimensional cones
using Hirzebruch-Jung continued fractions.

```@docs
cls_cone_normal_form
hirzebruch_jung
hilbert_basis
remove_vertex
```

## Computing subpolygons

Subpolygons can be either computed in memory or on disk using HDF5. The latter
is useful for computations where the amount of polygons will grow very large
(e.g. several billions), which would overload available memory. Note that even
in the latter case, some polygons will have to be kept in memory in order to do
equivalence checks with them later. In fact, to minimize memory usage further,
we will only keep the _hashes_ of those polygons and compare those. We use
`128`-bit hashes here, which makes the probability of a hash collision
negligible even for several billions of polygons.

```@docs
SubpolygonStorage
InMemorySubpolygonStoragePreferences
InMemorySubpolygonStorage
HDFSubpolygonStoragePreferences
HDFSubpolygonStorage
initialize_subpolygon_storage
subpolygons_single_step
subpolygons
restore_hdf_subpolygon_storage_status
```

