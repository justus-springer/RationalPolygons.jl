@doc raw"""
    SubpolygonStorage{T <: Integer}

An abstract supertype of storage options for computing subpolygons. There are
two subtypes `InMemorySubpolygonStorage` and `OnDiskSubpolygonStorage`. The
former keeps all subpolygons in memory, the latter delegeates their storage to
the disk.

"""
abstract type SubpolygonStorage{T <: Integer} end

rationality(st :: SubpolygonStorage{T}) where {T <: Integer} = st.rationality

number_of_interior_lattice_points(st :: SubpolygonStorage{T}) where {T <: Integer} = st.number_of_interior_lattice_points

total_count(st :: SubpolygonStorage{T}) where {T <: Integer} = st.total_count
