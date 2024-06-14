
function initialize_hdf_subpolygon_storage(
        f :: Union{HDF5.File, HDF5.Group},
        Ps :: Vector{<:RationalPolygon{T}};
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = false) where {T <: Integer}

    !isempty(Ps) || error("must provide a non-empty list of starting polygons")
    k = rationality(first(Ps))
    n = number_of_interior_lattice_points(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
    all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

    write_attribute(f, "primitive", primitive)
    write_attribute(f, "use_affine_normal_form", use_affine_normal_form)
    write_attribute(f, "rationality", k)
    write_attribute(f, "number_of_interior_lattice_points", n)
    write_attribute(f, "total_count", length(Ps))
    write_attribute(f, "last_completed_area", maximum(normalized_area.(Ps)) + 1)
    write_attribute(f, "is_finished", false)
    
    for P ∈ Ps
        a = normalized_area(P)
        N = number_of_vertices(P)
        if !haskey(f, "a$a")
            create_group(f, "a$a")
        end
        write_polygon_dataset(f, "a$a/n$N", [P])
    end

end


function subpolygons_single_step(
        f :: Union{HDF5.File, HDF5.Group};
        logging :: Bool = false,
        T :: Type{<:Integer} = Int)

    k = read_attribute(f, "rationality")
    primitive = read_attribute(f, "primitive")
    use_affine_normal_form = read_attribute(f, "use_affine_normal_form")
    last_completed_area = read_attribute(f, "last_completed_area")
    current_area = last_completed_area - 1

    while current_area >= 3
        haskey(f, "a$(current_area)") && break
        current_area -= 1
    end
    if current_area < 3
        attrs(f)["is_finished"] = true
        return 
    end
    current_area_group = f["a$(current_area)"]

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
            if !haskey(f, path)
                create_polygon_dataset(f, path, k, Qn; T)
            end
            old_Qs = Set(read_polygon_dataset(f, path))
            filter!(Q -> Q ∉ old_Qs, Qs)
            write_polygon_dataset(f, path, collect(Qs))
            attrs(f)["total_count"] += length(Qs)
        end

    end

    attrs(f)["last_completed_area"] = current_area
    total_count = read_attribute(f, "total_count")

    logging && @info "[Area = $current_area]. Running total: $total_count"

end

is_finished(f :: Union{HDF5.File, HDF5.Group}) = read_attribute(f, "is_finished")

function subpolygons(
        f :: Union{HDF5.File, HDF5.Group};
        logging :: Bool = false)

    while !is_finished(f)
        subpolygons_single_step(f; logging)
    end

    return f

end
