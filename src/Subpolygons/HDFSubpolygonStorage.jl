mutable struct HDFSubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    hdf_path :: String
    hdf_group :: String
    hash_sets :: Dict{T,Set{UInt64}}
    HDFSubpolygonStorage{T}(hdf_path :: String, hdf_group :: String = "/") where {T <: Integer} =
    new{T}(hdf_path, hdf_group, Dict{T,Set{UInt64}}())
end

function initialize_hdf_subpolygon_storage(
        st :: HDFSubpolygonStorage{T},
        Ps :: Vector{<:RationalPolygon{T}};
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = false) where {T <: Integer}


    !isempty(Ps) || error("must provide a non-empty list of starting polygons")
    k = rationality(first(Ps))
    n = number_of_interior_lattice_points(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
    all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

    f = h5open(st.hdf_path, "cw"; swmr = true)
    if !haskey(f, st.hdf_group)
        g = create_group(f, st.hdf_group)
    else
        g = f[st.hdf_group]
    end

    write_attribute(g, "primitive", primitive)
    write_attribute(g, "use_affine_normal_form", use_affine_normal_form)
    write_attribute(g, "rationality", k)
    write_attribute(g, "number_of_interior_lattice_points", n)
    write_attribute(g, "total_count", length(Ps))
    write_attribute(g, "last_completed_area", maximum(normalized_area.(Ps)) + 1)
    write_attribute(g, "is_finished", false)
    
    for a = 3 : maximum(normalized_area.(Ps))
        st.hash_sets[a] = Set{UInt64}()
    end

    for P ∈ Ps
        a = normalized_area(P)
        N = number_of_vertices(P)
        if !haskey(g, "a$a")
            create_group(g, "a$a")
        end
        write_polygon_dataset(g, "a$a/n$N", [P])
        push!(st.hash_sets[a], hash(P))
    end

    close(f)

end

function restore_hash_set(st :: HDFSubpolygonStorage{T}) where {T <: Integer}

    f = h5open(st.hdf_path, "r+"; swmr = true)
    g = f[st.hdf_group]
    k = read_attribute(g, "rationality")

    st.hash_sets = Dict{T,Set{UInt64}}()
    last_completed_area = read_attribute(g, "last_completed_area")
    for a = 3 : last_completed_area - 1
        st.hash_sets[a] = Set{UInt64}()
    end

    as = [parse(Int, s[2:end]) for s ∈ keys(g)]
    filter!(a -> a < last_completed_area, as)
    for a ∈ as
        st.hash_sets[a] = Set{UInt64}()
        for n ∈ keys(g["a$a"])
            Ps = read_polygon_dataset(k, g["a$a"], n)
            for P ∈ Ps
                push!(st.hash_sets[a], hash(P))
            end
        end
    end

end


function subpolygons_single_step(
        st :: HDFSubpolygonStorage{T};
        logging :: Bool = false,
        block_size :: Int = 10^6) where {T <: Integer}

    f = h5open(st.hdf_path, "r+"; swmr = true)
    g = f[st.hdf_group]

    k = read_attribute(g, "rationality")
    primitive = read_attribute(g, "primitive")
    use_affine_normal_form = read_attribute(g, "use_affine_normal_form")
    last_completed_area = read_attribute(g, "last_completed_area")
    current_area = last_completed_area - 1

    while current_area >= 3
        haskey(g, "a$(current_area)") && break
        current_area -= 1
    end
    if current_area < 3
        attrs(g)["is_finished"] = true
        return 
    end
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

            logging && @info "[a = $(current_area), n = $n, block $b/$number_of_blocks]. Polygons to peel: $(length(Ps))."

            Threads.@threads for P ∈ Ps
                id = Threads.threadid()
                for j = 1 : n
                    Q, keeps_genus = remove_vertex(P, j; primitive)
                    keeps_genus || continue
                    Qn = number_of_vertices(Q)
                    Qn > 2 || continue
                    if use_affine_normal_form
                        Q = affine_normal_form(Q)
                    else
                        Q = unimodular_normal_form(Q)
                    end

                    Qa = normalized_area(Q)
                    if !haskey(out_dicts[id], (Qa,Qn))
                        out_dicts[id][(Qa,Qn)] = Set{RationalPolygon{T,Qn,2*Qn}}()
                    end
                    Q ∉ out_dicts[id][(Qa,Qn)] || continue
                    hash(Q) ∉ st.hash_sets[Qa] || continue
                    push!(out_dicts[id][(Qa,Qn)], Q)
                end
            end

            new_polygons_count = sum([sum(length.(values(d))) for d ∈ out_dicts])
            logging && @info "[a = $(current_area), n = $n]. Peeling complete. New polygons: $new_polygons_count."

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
                attrs(g)["total_count"] += length(Qs)
            end

            total_count = read_attribute(g, "total_count")
            hash_set_size = sum([length(Qs) for (_,Qs) ∈ st.hash_sets])
            hash_set_bytes = Base.format_bytes(Base.summarysize(st.hash_sets))
            
            logging && @info "[a = $current_area, n = $n]. Writeout complete. Running total: $total_count. Hash set size: $hash_set_size ($hash_set_bytes)."

        end

    end

    attrs(g)["last_completed_area"] = current_area
    delete!(st.hash_sets, current_area)

    close(f)

end

function is_finished(st :: HDFSubpolygonStorage{T}) where {T <: Integer}
    f = h5open(st.hdf_path, "r"; swmr = true)
    is_finished = read_attribute(f[st.hdf_group], "is_finished")
    close(f)
    return is_finished
end

function subpolygons(
        st :: HDFSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    while !is_finished(st)
        subpolygons_single_step(st; logging)
    end

    return st

end
