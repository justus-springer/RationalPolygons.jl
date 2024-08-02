@doc raw"""
    SubpolygonStorage{T <: Integer}

An abstract supertype of storage options for computing subpolygons. There are
two subtypes `InMemorySubpolygonStorage` and `HDFSubpolygonStorage`. The
former keeps all subpolygons in memory, the latter delegeates their storage to
the disk using the HDF5 format. Both implement `subpolygons_single_step`.

"""
abstract type SubpolygonStorage{T <: Integer} end


@doc raw"""
    subpolygons(st :: SubpolygonStorage{T}; logging :: Bool = false) where {T <: Integer}

Compute all subpolygons with the given storage option.

"""
function subpolygons(st :: SubpolygonStorage{T}; logging :: Bool = false) where {T <: Integer}

    while last_completed_area(st) > 1
        subpolygons_single_step(st; logging)
    end

    return st

end
