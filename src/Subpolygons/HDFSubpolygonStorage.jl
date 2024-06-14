struct HDFSubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    hdf_path :: String
    hdf_group :: String
    HDFSubpolygonStorage{T}(hdf_path :: String, hdf_group :: String = "/") where {T <: Integer} = new{T}(hdf_path, hdf_group)
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

    f = h5open(st.hdf_path, "cw")
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
    
    for P ∈ Ps
        a = normalized_area(P)
        N = number_of_vertices(P)
        if !haskey(g, "a$a")
            create_group(g, "a$a")
        end
        write_polygon_dataset(g, "a$a/n$N", [P])
    end

end


function subpolygons_single_step(
        st :: HDFSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    f = h5open(st.hdf_path, "r+")
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

        Ps = read_polygon_dataset(current_area_group, "n$n") 

        out_dicts = Dict{Tuple{T,Int}}{Set{<:RationalPolygon{T}}}[]
        for i = 1 : Threads.nthreads()
            push!(out_dicts, Dict{Tuple{T,Int}}{Set{<:RationalPolygon{T}}}())
        end

        logging && @info "[Area = $(current_area), #vertices = $n]. Polygons to peel: $(length(Ps))."

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
                push!(out_dicts[id][(Qa,Qn)], Q)
            end
        end

        total_dict = mergewith(union!, out_dicts...)

        for ((Qa,Qn),Qs) ∈ total_dict
            path = "a$(Qa)/n$(Qn)"
            if !haskey(g, path)
                create_polygon_dataset(g, path, k, Qn; T)
            end
            old_Qs = Set(read_polygon_dataset(g, path))
            filter!(Q -> Q ∉ old_Qs, Qs)
            write_polygon_dataset(g, path, collect(Qs))
            attrs(g)["total_count"] += length(Qs)
        end

    end

    attrs(g)["last_completed_area"] = current_area
    total_count = read_attribute(g, "total_count")

    logging && @info "[Area = $current_area]. Running total: $total_count"

    close(f)

end

function is_finished(st :: HDFSubpolygonStorage{T}) where {T <: Integer}
    f = h5open(st.hdf_path, "r")
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
