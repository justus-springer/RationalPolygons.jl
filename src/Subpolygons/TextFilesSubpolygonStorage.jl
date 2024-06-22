struct TextFilesSubpolygonsPreferences{T <: Integer}
    rationality :: T
    number_of_interior_lattice_points :: Int
    primitive :: Bool
    use_affine_normal_form :: Bool

    TextFilesSubpolygonsPreferences{T}(;
        rationality :: T = one(T),
        number_of_interior_lattice_points :: Int = 1,
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = false) where {T <: Integer} =
    new{T}(rationality, number_of_interior_lattice_points, primitive, use_affine_normal_form)

end

mutable struct TextFilesSubpolygonStorage{T} <: SubpolygonStorage{T}
    preferences :: TextFilesSubpolygonsPreferences{T}
    directory :: String
    hash_sets :: Dict{T,Set{UInt128}}
    last_completed_area :: T
    total_count :: Int

    function TextFilesSubpolygonStorage{T}(preferences :: TextFilesSubpolygonsPreferences{T}, directory :: String) where {T <: Integer}
        isdir(directory) || error("$directory is not a directory")
        return new{T}(preferences, directory, Dict{T,Set{UInt128}}(), 0, 0)
    end

    TextFilesSubpolygonStorage{T}(directory :: String;
            rationality :: T = one(T),
            number_of_interior_lattice_points :: Int = 1,
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false) where {T <: Integer} = 
    TextFilesSubpolygonStorage{T}(TextFilesSubpolygonsPreferences{T}(;rationality, number_of_interior_lattice_points, primitive, use_affine_normal_form), directory)

    function TextFilesSubpolygonStorage{T}(
            Ps :: Vector{<:RationalPolygon{T}},
            directory :: String;
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false) where {T <: Integer}
        k = rationality(first(Ps))
        n = number_of_interior_lattice_points(first(Ps))
        all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

        pref = TextFilesSubpolygonsPreferences{T}(;rationality = k, number_of_interior_lattice_points = n, primitive, use_affine_normal_form)
        st = TextFilesSubpolygonStorage{T}(pref, directory)
        initialize_subpolygon_storage(st, Ps)

    end

end

function export_text_files_subpolygon_storage_status(st :: TextFilesSubpolygonStorage{T}) where {T <: Integer}
    open(joinpath(st.directory, "last_completed_area.txt"), "w") do f
        println(f, st.last_completed_area)
    end
end

function restore_text_files_subpolygon_storage_status(st :: TextFilesSubpolygonStorage{T}) where {T <: Integer}
    st.last_completed_area = parse(Int, first(readlines(joinpath(st.directory, "last_completed_area.txt"))))
    st.total_count = 0
    for filename ∈ readdir(st.directory)
        startswith(filename, "a") || continue
        filepath = joinpath(st.directory, filename)
        a = parse(Int, filename[2:end-4])
        st.total_count += countlines(filepath)
        a < st.last_completed_area || continue
        Ps = parse_rational_polygons(st.preferences.rationality, filepath)
        st.hash_sets[a] = Set([xxh3_128(vertex_matrix(P)) for P ∈ Ps])
    end

end

function initialize_subpolygon_storage(st :: TextFilesSubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    isempty(readdir(st.directory)) || error("$(st.directory) is dirty. Please provide an empty directory")
    maximum_area = maximum(normalized_area.(Ps))
    for a = 3 : maximum_area
        touch(joinpath(st.directory, "a$a.txt"))
        st.hash_sets[a] = Set{UInt128}()
    end
    for P ∈ Ps
        a = normalized_area(P)
        filepath = joinpath(st.directory, "a$a.txt")
        write_rational_polygons([P], filepath, "a")
        push!(st.hash_sets[a], xxh3_128(vertex_matrix(P)))
    end
    st.last_completed_area = maximum_area + 1
    st.total_count = length(Ps)
    export_text_files_subpolygon_storage_status(st)

    return st

end

is_finished(st :: TextFilesSubpolygonStorage{T}) where {T <: Integer} = st.last_completed_area <= 3

function subpolygons_single_step(
        st :: TextFilesSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    current_area = st.last_completed_area - 1
    filepath = joinpath(st.directory, "a$current_area.txt")
    Ps = parse_rational_polygons(st.preferences.rationality, filepath)

    logging && @info "[a = $current_area]. Polygons to peel: $(length(Ps))."

    out_array = Dict{T, Set{RationalPolygon{T}}}[]
    for i = 1 : Threads.nthreads()
        push!(out_array, Dict{T, Set{RationalPolygon{T}}}())
    end

    Threads.@threads for P ∈ Ps
        id = Threads.threadid()
        for j = 1 : number_of_vertices(P)
            Q, keeps_genus = remove_vertex(P, j; st.preferences.primitive)
            keeps_genus || continue
            number_of_vertices(Q) > 2 || continue
            Q = st.preferences.use_affine_normal_form ? affine_normal_form(Q) : unimodular_normal_form(Q)
            a = normalized_area(Q)
            if !haskey(out_array[id], a)
                out_array[id][a] = Set{RationalPolygon{T}}()
            end
            xxh3_128(vertex_matrix(Q)) ∉ st.hash_sets[a] || continue
            push!(out_array[id][a], Q)
        end
    end

    new_polygons = mergewith(union!, out_array...)
    new_polygons_count = sum(length.(values(new_polygons)))
    logging && @info "[a = $current_area]. Peeling complete. New polygons: $new_polygons_count"

    for (a,Qs) ∈ new_polygons
        filepath = joinpath(st.directory, "a$a.txt")
        write_rational_polygons(collect(Qs), filepath, "a")
        union!(st.hash_sets[a], xxh3_128.(Qs))
        st.total_count += length(Qs)
    end

    st.last_completed_area = current_area
    delete!(st.hash_sets, current_area)
    export_text_files_subpolygon_storage_status(st)

    logging && @info "[a = $current_area]. Writeout complete. Running total: $(st.total_count)"

end

function subpolygons(
        st :: TextFilesSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    while !is_finished(st)
        subpolygons_single_step(st; logging)
    end

    return st

end
