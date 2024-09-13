@doc raw"""
    struct HDFSubpolygonStoragePreferences{T <: Integer}

A struct holding preferences for `HDFSubpolygonStorage`. There are the
following fields:

- `rationality :: T`: The rationality of the polygons. The default is `one(T)`.
- `primitive :: Bool`: Whether only subpolygons with primitive vertices should
   be computed. The default is `false`.
- `use_affine_normal_form :: Bool`: Whether to use [`affine_normal_form`](@ref)
    or [`unimodular_normal_form`](@ref). The default is `true`, i.e. affine normal
    form.
- `only_equal_number_of_interior_lattice_points :: Bool`: Whether only
    subpolygons having the same number of interior lattice points as the starting
    polygons should be computed. The default is `false`.
- `exclude_very_thin_polygons`: Whether polygons that can be realized in ``\mathbb{R} \times [0,1]`` should be excluded. This is only relevant for polygons with no interior lattice points. The default is `false`.
- `block_size :: Int`: How many polygons should be read into memory at once
    during the shaving process. Defaults to `10^6`.
- `maximum_number_of_vertices :: Int`: An upper bound for the maximal number of
    vertices to be expected in the computation. This has to be set since every
    HDF5 file generated will have a dataset "numbers\_of\_polygons" storing the
    number of polygons for each number of vertices and the size of this dataset
    needs to be set beforehand. Defaults to `100`, which should be more than enough for any
    feasable computation.
- `swmr :: Bool`: Whether to use single-reader-multiple-writer mode for HDF5.
    Defaults to `true`.

"""
struct HDFSubpolygonStoragePreferences{T <: Integer}
    rationality :: T
    primitive :: Bool
    use_affine_normal_form :: Bool
    only_equal_number_of_interior_lattice_points :: Bool
    exclude_very_thin_polygons :: Bool
    block_size :: Int
    swmr :: Bool
    maximum_number_of_vertices :: Int

    HDFSubpolygonStoragePreferences{T}(;
        rationality :: T = one(T),
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = true,
        only_equal_number_of_interior_lattice_points :: Bool = false,
        exclude_very_thin_polygons :: Bool = false,
        block_size :: Int = 10^6,
        swmr :: Bool = true,
        maximum_number_of_vertices :: Int = 100) where {T <: Integer} =
    new{T}(rationality, primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons, block_size, swmr, maximum_number_of_vertices)

end


@doc raw"""
    mutable struct HDFSubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}   

A struct for managing the result of a subpolygon computation using the HDF5
file format. It has the following fields:

- `preferences :: HDFSubpolygonStoragePreferences{T}`
- `file_path :: String`: The file path of the HDF5 file to be created.
- `group_path :: String`: Path to a group within the HDF5 file, if it already
   exists. Defaults to "/", i.e. the root group.
- `hash_sets :: Dict{T,Set{UInt128}}`: A dictionary of hashes of polygons
   already encountered. This is used for comparison with new polygons to ensure
   the result contains each polygon exactly once. We hold these hashes in memory
   to avoid needing to read in polygons that have been written out in the past,
   which saves a lot of time.
- `last_completed_area :: T`: The last area that has been completed. This
   counts down from the maximum area of the starting polygons to 1.
- `total_count :: Int`

"""
mutable struct HDFSubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    preferences :: HDFSubpolygonStoragePreferences{T}
    file_path :: String
    group_path :: String
    hash_sets :: Dict{T,Set{UInt128}}
    last_completed_area :: T
    total_count :: Int


    HDFSubpolygonStorage{T}(preferences :: HDFSubpolygonStoragePreferences{T}, file_path :: String, group_path :: String = "/") where {T <: Integer} =
    new{T}(preferences, file_path, group_path, Dict{T,Set{UInt128}}(), 0, 0)

    HDFSubpolygonStorage{T}(file_path :: String, 
            group_path :: String = "/";
            rationality :: T = one(T),
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = true,
            only_equal_number_of_interior_lattice_points :: Bool = false,
            exclude_very_thin_polygons :: Bool = false,
            block_size :: Int = 10^6,
            swmr :: Bool = true,
            maximum_number_of_vertices :: Int = 100) where {T <: Integer} = 
    HDFSubpolygonStorage{T}(HDFSubpolygonStoragePreferences{T}(;rationality, primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons, block_size, swmr, maximum_number_of_vertices), file_path, group_path)

    function HDFSubpolygonStorage{T}(Ps :: Vector{<:RationalPolygon{T}},
            file_path :: String,
            group_path :: String = "/";
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false,
            only_equal_number_of_interior_lattice_points :: Bool = true,
            exclude_very_thin_polygons :: Bool = false,
            block_size :: Int = 10^6,
            swmr :: Bool = true,
            maximum_number_of_vertices :: Int = 100) where {T <: Integer}

        k = rationality(first(Ps))
        all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")

        pref = HDFSubpolygonStoragePreferences{T}(;rationality = k, primitive, use_affine_normal_form, only_equal_number_of_interior_lattice_points, exclude_very_thin_polygons, block_size, swmr, maximum_number_of_vertices)
        st = HDFSubpolygonStorage{T}(pref, file_path, group_path)
        initialize_subpolygon_storage(st, Ps)

    end

end

last_completed_area(st :: HDFSubpolygonStorage) = st.last_completed_area

@doc raw"""
    initialize_subpolygon_storage(st :: HDFSubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

Initialize a newly created `HDFSubpolygonStorage` with a given set of starting
polygons. To restart an interrupted computation, see also
[`restore_hdf_subpolygon_storage_status`](@ref).

"""
function initialize_subpolygon_storage(st :: HDFSubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    k = rationality(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")

    f = h5open(st.file_path, "cw"; swmr = st.preferences.swmr)
    if !haskey(f, st.group_path)
        g = create_group(f, st.group_path)
    else
        g = f[st.group_path]
    end

    maximum_area = maximum(normalized_area.(Ps))
    st.last_completed_area = maximum_area + 1
    attrs(g)["last_completed_area"] = st.last_completed_area
    
    for a = 1 : maximum_area
        st.hash_sets[a] = Set{UInt128}()
        create_group(g, "a$a")
    end

    st.total_count = length(Ps)
    for P ∈ Ps
        a = normalized_area(P)
        N = number_of_vertices(P)
        write_polygon_dataset(g, "a$a/n$N", [P])
        push!(st.hash_sets[a], xxh3_128(vertex_matrix(P)))
    end

    g["numbers_of_polygons"] = zeros(Int, maximum_area, st.preferences.maximum_number_of_vertices) 
    numbers_of_polygons = g["numbers_of_polygons"]
    for a = 1 : maximum_area
        for n = 1 : st.preferences.maximum_number_of_vertices
            haskey(g["a$a"], "n$n") || continue
            val = length(dataspace(g["a$a"]["n$n"]))
            numbers_of_polygons[a,n] = val
        end
    end

    close(f)

    return st

end


@doc raw"""
    restore_hdf_subpolygon_storage_status(st :: HDFSubpolygonStorage{T}) where {T <: Integer}

Restore a subpolygons computation's status that was interrupted from an HDF5
file. All polygons will be read in, hashed and saved into `st.hash_sets`.
Moreover, `st.last_completed_area` and `st.total_count` will be properly set.
After calling this function, the computation can be resumed by calling
`subpolygons(st)` again.

"""
function restore_hdf_subpolygon_storage_status(st :: HDFSubpolygonStorage{T}) where {T <: Integer}

    f = h5open(st.file_path, "r+"; swmr = st.preferences.swmr)
    g = f[st.group_path]
    k = st.preferences.rationality

    st.last_completed_area = read_attribute(g, "last_completed_area")

    A = read_dataset(g, "numbers_of_polygons")
    st.total_count = sum(A)
    for a_key ∈ keys(g)
        startswith(a_key, "a") || continue
        a = parse(Int, a_key[2:end])
        st.hash_sets[a] = Set{UInt128}()
        for n_key ∈ keys(g[a_key])
            startswith(n_key, "n") || continue
            n = parse(Int, n_key[2:end])
            HDF5.set_extent_dims(g[a_key][n_key], (A[a,n],))

            Ps = read_polygon_dataset(k, g[a_key], n_key)
            a < st.last_completed_area || continue
            for P ∈ Ps
                push!(st.hash_sets[a], xxh3_128(vertex_matrix(P)))
            end
        end
    end

    close(f)

end


@doc raw"""
    subpolygons_single_step(st :: HDFSubpolygonStorage{T}; logging :: Bool = false) where {T <: Integer}

Perform a single step in a subpolygon computation using the HDF5 file format.

"""
function subpolygons_single_step(st :: HDFSubpolygonStorage{T}; logging :: Bool = false) where {T <: Integer}

    f = h5open(st.file_path, "r+"; swmr = st.preferences.swmr)
    g = f[st.group_path]
    k = st.preferences.rationality
    block_size = st.preferences.block_size

    current_area = read_attribute(g, "last_completed_area") - 1
    current_area_group = g["a$(current_area)"]

    ns = sort([parse(Int, s[2:end]) for s ∈ keys(current_area_group)])
    for n ∈ ns

        N = length(dataspace(current_area_group["n$n"]))
        number_of_blocks = N ÷ block_size + 1

        for b = 1 : number_of_blocks

            I = (b-1)*block_size + 1 : min(b*block_size, N)

            Ps = read_polygon_dataset(k, current_area_group, "n$n", I) 

            out_dicts = Dict{Tuple{T,Int}}{Set{<:RationalPolygon{T}}}[]
            for i = 1 : Threads.nthreads()
                push!(out_dicts, Dict{Tuple{T,Int}}{Set{<:RationalPolygon{T}}}())
            end

            logging && @info "[a = $current_area, n = $n, block $b/$number_of_blocks]. Polygons to peel: $(length(Ps))."

            Threads.@threads for P ∈ Ps
                id = Threads.threadid()
                for j = 1 : n
                    Q, keeps_genus = remove_vertex(P, j; primitive = st.preferences.primitive)
                    if st.preferences.only_equal_number_of_interior_lattice_points 
                        keeps_genus || continue
                    end
                    Qn = number_of_vertices(Q)
                    Qn > 2 || continue

                    if st.preferences.exclude_very_thin_polygons
                        minimal_number_of_interior_integral_lines(Q) > 0 || continue
                    end

                    Q = st.preferences.use_affine_normal_form ? affine_normal_form(Q) : unimodular_normal_form(Q)
                    Qa = normalized_area(Q)
                    if !haskey(out_dicts[id], (Qa,Qn))
                        out_dicts[id][(Qa,Qn)] = Set{RationalPolygon{T,Qn,2*Qn}}()
                    end
                    xxh3_128(vertex_matrix(Q)) ∉ st.hash_sets[Qa] || continue
                    push!(out_dicts[id][(Qa,Qn)], Q)
                end
            end

            new_polygons_count = sum([sum(length.(values(d))) for d ∈ out_dicts])
            logging && @info "[a = $(current_area), n = $n, block $b/$number_of_blocks]. Peeling complete. New polygons: $new_polygons_count."

            for i = 2 : length(out_dicts)
                mergewith!(union!, out_dicts[1], out_dicts[i])
                out_dicts[i] = Dict{Tuple{T,Int}}{Set{<:RationalPolygon{T}}}()
            end

            for ((Qa,Qn), Qs) ∈ out_dicts[1]
                path = "a$(Qa)/n$(Qn)"
                write_polygon_dataset(g, path, collect(Qs))
                for Q ∈ Qs
                    push!(st.hash_sets[Qa], xxh3_128(vertex_matrix(Q)))
                end
                g["numbers_of_polygons"][Qa,Qn] += length(Qs)
                st.total_count += length(Qs)
            end

            hash_set_size = sum([length(Qs) for (_,Qs) ∈ st.hash_sets])
            hash_set_bytes = Base.format_bytes(Base.summarysize(st.hash_sets))
            
            logging && @info "[a = $current_area, n = $n, block $b/$number_of_blocks]. Writeout complete. Running total: $(st.total_count). Hash set size: $hash_set_size ($hash_set_bytes)."

        end

    end

    st.last_completed_area = current_area
    delete!(st.hash_sets, current_area)
    attrs(g)["last_completed_area"] = current_area

    close(f)

end
