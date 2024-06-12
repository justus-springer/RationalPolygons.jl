
struct HDFSubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    hdf_file_path :: String
    group_path :: String
    HDFSubpolygonStorage{T}(hdf_file_path :: String, group_path :: String = "/") where {T <: Integer} =
    new{T}(hdf_file_path, group_path)
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

    f = h5open(st.hdf_file_path, "cw")
    if haskey(f, st.group_path)
        g = f[st.group_path]
    else
        g = create_group(f, st.group_path)
    end

    attrs(g)["primitive"] = primitive
    attrs(g)["use_affine_normal_form"] = use_affine_normal_form
    attrs(g)["rationality"] = k
    attrs(g)["number_of_interior_lattice_points"] = n
    
    for P ∈ Ps
        a = normalized_area(P)
        N = number_of_vertices(P)
        if !haskey(g, "a$a")
            create_group(g, "a$a")
        end
        if !haskey(g["a$a"], "n$N")
            export_polygons_to_h5(g["a$a"], "n$N", [P])
        else
            append_polygons_to_h5(g["a$a"], "n$N", [P])
        end
    end

    attrs(g)["total_count"] = length(Ps)
    attrs(g)["last_completed_area"] = maximum(normalized_area.(Ps)) + 1
    attrs(g)["is_finished"] = false

    close(f)

end


function subpolygons_single_step(
        st :: HDFSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    f = h5open(st.hdf_file_path, "r+")
    g = f[st.group_path]
    k = attrs(g)["rationality"]
    primitive = attrs(g)["primitive"]
    use_affine_normal_form = attrs(g)["use_affine_normal_form"]

    last_completed_area = attrs(g)["last_completed_area"]
    current_area = last_completed_area - 1
    while current_area >= 3
        haskey(g, "a$(current_area)") && break
        current_area -= 1
    end
    if current_area < 3
        attrs(g)["is_finished"] = true
        total_count = attrs(g)["total_count"]
        close(f)
        return total_count
    end
    current_area_group = g["a$(current_area)"]


    for n_string ∈ keys(current_area_group)

        n = parse(Int,n_string[2:end])
        Ps = import_polygons_from_h5(k, current_area_group, n_string) 

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

        for ((Qa,Qn), Qs) ∈ total_dict
            if !haskey(g, "a$(Qa)/n$(Qn)")
                export_polygons_to_h5(g, "a$(Qa)/n$(Qn)", collect(Qs))
                attrs(g)["total_count"] += length(Qs)
            else
                old_Qs = Set(import_polygons_from_h5(k, g, "a$(Qa)/n$(Qn)"))
                filter!(Q -> Q ∉ old_Qs, Qs)
                append_polygons_to_h5(g, "a$(Qa)/n$(Qn)", collect(Qs))
                attrs(g)["total_count"] += length(Qs)
            end
        end

    end

    attrs(g)["last_completed_area"] = current_area
    total_count = attrs(g)["total_count"]
    close(f)

    logging && @info "[Area = $current_area]. Running total: $total_count"

    return total_count

end

function is_finished(st :: HDFSubpolygonStorage{T}) where {T <: Integer}
    f = h5open(st.hdf_file_path, "r+")
    g = f[st.group_path]
    is_finished = attrs(g)["is_finished"]
    close(f)
    return is_finished
end

return_value(::HDFSubpolygonStorage) = nothing
