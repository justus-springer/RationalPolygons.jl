@doc raw"""
    struct InMemorySubpolygonStoragePreferences{T <: Integer}

A struct holding preferences for `InMemorySubpolygonStorage`. There are the following fields:

- `primitive :: Bool`: Whether only subpolygons with primitive vertices should
   be computed. The default is `false`.
- `use_affine_normal_form :: Bool`: Whether to use [`affine_normal_form`](@ref)
    or [`unimodular_normal_form`](@ref). The default is `true`, i.e. affine normal
    form.
- `only_equal_number_of_interior_lattice_points :: Bool`: Whether only
    subpolygons having the same number of interior lattice points as the starting
    polygons should be computed. The default is `false`.
- `exclude_very_thin_polygons`: Whether polygons that can be realized in ``\mathbb{R} \times [0,1]`` should be excluded. This is only relevant for polygons with no interior lattice points. The default is `false`.

"""
struct InMemorySubpolygonStoragePreferences{T <: Integer} 
    primitive :: Bool
    use_affine_normal_form :: Bool
    only_equal_number_of_interior_lattice_points :: Bool
    exclude_very_thin_polygons :: Bool

    InMemorySubpolygonStoragePreferences{T}(;
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = true,
        only_equal_number_of_interior_lattice_points :: Bool = false,
        exclude_very_thin_polygons :: Bool = false) where {T <: Integer} =
    new{T}(primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons)

end


@doc raw"""
    mutable struct InMemorySubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}   

A struct holding results of a subpolygon computation. It has the following fields:

- `preferences :: InMemorySubpolygonStoragePreferences{T}`
- `polygons :: Dict{T,Set{RationalPolygon{T}}}`
- `last_completed_area :: T`
- `total_count :: Int`

"""
mutable struct InMemorySubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    preferences :: InMemorySubpolygonStoragePreferences{T}
    polygons :: Dict{T,Set{RationalPolygon{T}}}
    last_completed_area :: T
    total_count :: Int

    InMemorySubpolygonStorage{T}(preferences :: InMemorySubpolygonStoragePreferences{T}) where {T <: Integer} =
    new{T}(preferences, Dict{T,Set{RationalPolygon{T}}}(), 0, 0)

    InMemorySubpolygonStorage{T}(;
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = true,
            only_equal_number_of_interior_lattice_points  :: Bool = false,
            exclude_very_thin_polygons :: Bool = false) where {T <: Integer} =
    InMemorySubpolygonStorage{T}(InMemorySubpolygonStoragePreferences{T}(;primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons))

    function InMemorySubpolygonStorage{T}(
            Ps :: Vector{<:RationalPolygon{T}};
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = true,
            only_equal_number_of_interior_lattice_points :: Bool = false,
            exclude_very_thin_polygons :: Bool = false) where {T <: Integer}

        pref = InMemorySubpolygonStoragePreferences{T}(;primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons)
        st = InMemorySubpolygonStorage{T}(pref)
        initialize_subpolygon_storage(st, Ps)

    end

end

last_completed_area(st :: InMemorySubpolygonStorage) = st.last_completed_area

function initialize_subpolygon_storage(st :: InMemorySubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    k = rationality(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")

    maximum_area = maximum(normalized_area.(Ps))
    for a = 1 : maximum_area
        st.polygons[a] = Set{RationalPolygon{T}}()
    end
    for P ∈ Ps
        push!(st.polygons[normalized_area(P)], P)
    end
    st.last_completed_area = maximum_area + 1
    st.total_count = length(Ps)

    return st

end


@doc raw"""
    subpolygons_single_step(st :: InMemorySubpolygonStorage{T}; logging :: Bool = false) where {T <: Integer}

Perform a single step of a subpolygon computation in memory.

"""
function subpolygons_single_step(st :: InMemorySubpolygonStorage{T}; logging :: Bool = false) where {T <: Integer}

    current_area = st.last_completed_area - 1
    Ps = collect(st.polygons[current_area])

    logging && @info "[a = $current_area]. Polygons to peel: $(length(Ps))."

    out_array = Set{RationalPolygon{T}}[]
    for i = 1 : Threads.nthreads()
        push!(out_array, Set{RationalPolygon{T}}())
    end

    Threads.@threads for P ∈ Ps
        for j = 1 : number_of_vertices(P)
            Q, keeps_genus = remove_vertex(P, j; primitive =  st.preferences.primitive)
            if st.preferences.only_equal_number_of_interior_lattice_points 
                keeps_genus || continue
            end

            number_of_vertices(Q) > 2 || continue

            if st.preferences.exclude_very_thin_polygons
                minimal_number_of_interior_integral_lines(Q) > 0 || continue
            end

            if st.preferences.use_affine_normal_form
                Q = affine_normal_form(Q)
            else
                Q = unimodular_normal_form(Q)
            end
            Q ∉ st.polygons[normalized_area(Q)] || continue
            push!(out_array[Threads.threadid()], Q)
        end
    end

    new_polygons = union!(out_array...)

    for P ∈ new_polygons
        a = normalized_area(P)
        push!(st.polygons[a], P)
        st.total_count += 1
    end

    st.last_completed_area = current_area

    logging && @info "[a = $current_area]. Peeling complete. New polygons: $(length(new_polygons)). Running total: $(st.total_count)"

end


@doc raw"""
    subpolygons(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    subpolygons(P :: RationalPolygon{T}) where {T <: Integer}

Compute all subpolygons of a rational polygon or list of rational polygons. The
computation is done in memory, for storage on disk see also
[`HDFSubpolygonStorage`](@ref). This function takes the following keyword
arguments:

- `primitive :: Bool`: Whether only subpolygons with primitive vertices should
    be returned. The default is `false`.
- `use_affine_normal_form :: Bool`: Whether to use affine or unimodular normal
    form. The default is `true`, so affine normal form.
- `only_equal_number_of_interior_lattice_points :: Bool`: Whether only
    subpolygons that share the same number of interior lattice points with the
    starting polygons should be returned.
- `logging :: Bool`: Whether to display logging messages about the computation
    progress.

# Example

There are 148 subpolygons of the square of side length 3, up to affine
equivalence. The maximal number of vertices of those is 8, attained by exactly
one polygon.

```jldoctest
julia> Ps = subpolygons(convex_hull(LatticePoint{Int}[(0,0),(3,0),(0,3),(3,3)]));

julia> length(Ps)
148

julia> maximum(number_of_vertices.(Ps))
8
```

"""
function subpolygons(Ps :: Vector{<:RationalPolygon{T}}; 
    primitive :: Bool = false, 
    use_affine_normal_form :: Bool = true,
    only_equal_number_of_interior_lattice_points :: Bool = false,
    exclude_very_thin_polygons :: Bool = false,
    logging :: Bool = false) where {T <: Integer}

    st = InMemorySubpolygonStorage{T}(Ps; primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons)
    subpolygons(st; logging)
    return collect(union(values(st.polygons)...))
end

subpolygons(P :: RationalPolygon{T}; 
    primitive :: Bool = false, 
    use_affine_normal_form :: Bool = true,
    only_equal_number_of_interior_lattice_points :: Bool = false,
    exclude_very_thin_polygons :: Bool = false,
    logging :: Bool = false) where {T <: Integer} =
subpolygons([P]; primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons, logging)
