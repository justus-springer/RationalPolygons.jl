
mutable struct InMemorySubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    preferences :: SubpolygonsPreferences{T}
    polygons_dict :: Dict{T,Set{RationalPolygon{T}}}
    last_completed_area :: T
    total_count :: Int

    InMemorySubpolygonStorage{T}(preferences :: SubpolygonsPreferences{T}) where {T <: Integer} =
    new{T}(preferences, Dict{T,Set{RationalPolygon{T}}}(), 0, 0)

    InMemorySubpolygonStorage{T}(;
            rationality :: T = one(T),
            number_of_interior_lattice_points :: Int = 1,
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false) where {T <: Integer} = 
    InMemorySubpolygonStorage{T}(SubpolygonsPreferences{T}(rationality, number_of_interior_lattice_points, primitive, use_affine_normal_form))

    function InMemorySubpolygonStorage{T}(
            Ps :: Vector{<:RationalPolygon{T}};
            primitive :: Bool = false,
            use_affine_normal_form :: Bool = false) where {T <: Integer}
        k = rationality(first(Ps))
        n = number_of_interior_lattice_points(first(Ps))
        all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

        pref = SubpolygonsPreferences{T}(k, n, primitive, use_affine_normal_form)
        st = InMemorySubpolygonStorage{T}(pref)
        initialize_subpolygon_storage(st, Ps)

    end

end

function initialize_subpolygon_storage(st :: InMemorySubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    maximum_area = maximum(normalized_area.(Ps))
    for a = 3 : maximum_area
        st.polygons_dict[a] = Set{RationalPolygon{T}}()
    end
    for P ∈ Ps
        push!(st.polygons_dict[normalized_area(P)], P)
    end
    st.last_completed_area = maximum_area + 1
    st.total_count = length(Ps)

    return st

end

is_finished(st :: InMemorySubpolygonStorage{T}) where {T <: Integer} = st.last_completed_area <= 3


function subpolygons_single_step(
        st :: InMemorySubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    current_area = maximum(filter(b -> b < st.last_completed_area, keys(st.polygons_dict)))
    Ps = collect(st.polygons_dict[current_area])

    logging && @info "[a = $current_area]. Polygons to peel: $(length(Ps))."

    out_array = Set{RationalPolygon{T}}[]
    for i = 1 : Threads.nthreads()
        push!(out_array, Set{RationalPolygon{T}}())
    end

    Threads.@threads for P ∈ Ps
        for j = 1 : number_of_vertices(P)
            Q, keeps_genus = remove_vertex(P, j; st.preferences.primitive)
            keeps_genus || continue
            number_of_vertices(Q) > 2 || continue
            if st.preferences.use_affine_normal_form
                Q = affine_normal_form(Q)
            else
                Q = unimodular_normal_form(Q)
            end
            Q ∉ st.polygons_dict[normalized_area(Q)] || continue
            push!(out_array[Threads.threadid()], Q)
        end
    end

    new_polygons = union!(out_array...)
    logging && @info "[a = $current_area]. Peeling complete. New polygons: $(length(new_polygons))"

    for P ∈ new_polygons
        a = normalized_area(P)
        push!(st.polygons_dict[a], P)
        st.total_count += 1
    end

    st.last_completed_area = current_area

    logging && @info "[a = $current_area]. Writeout complete. Running total: $(st.total_count)"

end

function subpolygons(
        st :: InMemorySubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    while !is_finished(st)
        subpolygons_single_step(st; logging)
    end

    return collect(union(values(st.polygons_dict)...))

end

subpolygons(Ps :: Vector{<:RationalPolygon{T}}; 
    primitive :: Bool = false, 
    use_affine_normal_form :: Bool = false,
    logging :: Bool = false) where {T <: Integer} =
subpolygons(InMemorySubpolygonStorage{T}(Ps; primitive, use_affine_normal_form); logging)

subpolygons(P :: RationalPolygon{T}; 
    primitive :: Bool = false, 
    use_affine_normal_form :: Bool = false,
    logging :: Bool = false) where {T <: Integer} =
subpolygons([P]; primitive, use_affine_normal_form, logging)
