# Subpolygons

Given a ``k``-rational polygon ``P``, we want to find all subpolygons of ``P``
up to equivalence. The algorithm used by RationalPoylgons.jl is described in
section 2.3 of [BS24](@cite). The main idea is to succesively remove vertices
of ``P`` by computing hilbert bases.

## Hilbert bases

We follow [CLS11](@cite) to compute hilbert bases of two-dimensional cones
using Hirzebruch-Jung continued fractions.

```@docs
cls_cone_normal_form
hirzebruch_jung
hilbert_basis
remove_vertex
```

## Computing subpolygons

Subpolygons can be either computed in memory or on disk using HDF5. The
latter is useful for large computations, since the amount of data can easily
overload memory.

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

