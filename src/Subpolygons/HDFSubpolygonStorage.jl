struct HDFSubpolygonsPreferences{T <: Integer}
    rationality :: T
    number_of_interior_lattice_points :: Int
    primitive :: Bool
    use_affine_normal_form :: Bool
    block_size :: Int

    HDFSubpolygonsPreferences{T}(;
        rationality :: T = one(T),
        number_of_interior_lattice_points :: Int = 1,
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = false,
        block_size :: Int = 10^6) where {T <: Integer} =
    new{T}(rationality, number_of_interior_lattice_points, primitive, use_affine_normal_form, block_size)

end

mutable struct HDFSubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    preferences :: HDFSubpolygonsPreferences{T}
    file_path :: String
    group_path :: String
    hash_sets :: Dict{T,Set{UInt64}}
    last_completed_area :: T
    total_count :: Int


    HDFSubpolygonStorage{T}(preferences :: HDFSubpolygonsPreferences{T}, file_path :: String, group_path :: String = "/") where {T <: Integer} =
    new{T}(preferences, file_path, group_path, Dict{T,Set{UInt64}}(), 0, 0)

    HDFSubpolygonStorage{T}(file_path :: String, 
            group_path :: String = "/";
            rationality :: T = one(T),
            number_of_interior_lattice_points :: Int = 1,
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false,
            block_size :: Int = 10^6) where {T <: Integer} = 
    HDFSubpolygonStorage{T}(HDFSubpolygonsPreferences{T}(;rationality, number_of_interior_lattice_points, primitive, use_affine_normal_form, block_size), file_path, group_path)

    function HDFSubpolygonStorage{T}(Ps :: Vector{<:RationalPolygon{T}},
            file_path :: String,
            group_path :: String = "/";
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false,
            block_size :: Int = 10^6) where {T <: Integer}
        k = rationality(first(Ps))
        n = number_of_interior_lattice_points(first(Ps))
        all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

        pref = HDFSubpolygonsPreferences{T}(;rationality = k, number_of_interior_lattice_points = n, primitive, use_affine_normal_form, block_size)
        st = HDFSubpolygonStorage{T}(pref, file_path, group_path)
        initialize_subpolygon_storage(st, Ps)

    end

end

function initialize_subpolygon_storage(
        st :: HDFSubpolygonStorage{T},
        Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    k = rationality(first(Ps))
    n = number_of_interior_lattice_points(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
    all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

    f = h5open(st.file_path, "cw"; swmr = true)
    if !haskey(f, st.group_path)
        g = create_group(f, st.group_path)
    else
        g = f[st.group_path]
    end

    st.last_completed_area = maximum(normalized_area.(Ps)) + 1
    attrs(g)["last_completed_area"] = st.last_completed_area
    
    for a = 3 : st.last_completed_area - 1
        st.hash_sets[a] = Set{UInt64}()
        create_group(g, "a$a")
    end

    st.total_count = length(Ps)
    for P ∈ Ps
        a = normalized_area(P)
        N = number_of_vertices(P)
        write_polygon_dataset(g, "a$a/n$N", [P])
        push!(st.hash_sets[a], hash(P))
    end

    close(f)

    return st

end

function restore_hdf_subpolygon_storage_status(st :: HDFSubpolygonStorage{T}) where {T <: Integer}

    f = h5open(st.file_path, "r+"; swmr = true)
    g = f[st.group_path]
    k = st.preferences.rationality

    st.last_completed_area = read_attribute(g, "last_completed_area")

    for a_key ∈ keys(g)
        a = parse(Int, a_key[2:end])
        st.hash_sets[a] = Set{UInt64}()
        for n ∈ keys(g[a_key])
            st.total_count += length(dataspace(g[a_key][n]))
            Ps = read_polygon_dataset(k, g[a_key], n)
            a < st.last_completed_area || continue
            for P ∈ Ps
                push!(st.hash_sets[a], hash(P))
            end
        end
    end

    close(f)

end


function subpolygons_single_step(
        st :: HDFSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    f = h5open(st.file_path, "r+"; swmr = true)
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
                    Q, keeps_genus = remove_vertex(P, j; st.preferences.primitive)
                    keeps_genus || continue
                    Qn = number_of_vertices(Q)
                    Qn > 2 || continue
                    Q = st.preferences.use_affine_normal_form ? affine_normal_form(Q) : unimodular_normal_form(Q)
                    Qa = normalized_area(Q)
                    if !haskey(out_dicts[id], (Qa,Qn))
                        out_dicts[id][(Qa,Qn)] = Set{RationalPolygon{T,Qn,2*Qn}}()
                    end
                    hash(Q) ∉ st.hash_sets[Qa] || continue
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
                    push!(st.hash_sets[Qa], hash(Q))
                end
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

is_finished(st :: HDFSubpolygonStorage{T}) where {T <: Integer} = st.last_completed_area <= 3

function subpolygons(
        st :: HDFSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    while !is_finished(st)
        subpolygons_single_step(st; logging)
    end

    return st

end
