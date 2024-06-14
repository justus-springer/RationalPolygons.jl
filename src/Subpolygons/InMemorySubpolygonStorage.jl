
mutable struct InMemorySubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    rationality :: T
    number_of_interior_lattice_points :: Int
    primitive :: Bool
    use_affine_normal_form :: Bool
    polygons_dict :: Dict{T,Set{RationalPolygon{T}}}
    last_volume :: T
    total_count :: Int

    @doc raw"""
        InMemorySubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    Construct an in-memory storage container of subpolygons starting with the
    given list of polygons. When using this constructor, the subpolygons are
    not yet computed. To compute the subpolygons, use [`subpolygons`](@ref).
    
    """
    function InMemorySubpolygonStorage{T}(
       starting_polygons :: Vector{<:RationalPolygon{T}};
       primitive :: Bool = false,
       use_affine_normal_form :: Bool = false) where {T <: Integer}

        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")

        k = rationality(first(starting_polygons))
        n = number_of_interior_lattice_points(first(starting_polygons))
        all(P -> rationality(P) == k, starting_polygons) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, starting_polygons) || error("all polygons must have the same number of interior lattice points")


        polygons_dict = Dict{T, Set{RationalPolygon{T}}}()
        last_volume = maximum(normalized_area.(starting_polygons)) + 1
        total_count = 0 
        st = new{T}(k, n, primitive, use_affine_normal_form, polygons_dict, last_volume, total_count)
        save!(st, starting_polygons)
        return st
    end

end

function save!(st :: InMemorySubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    for P ∈ Ps
        a = normalized_area(P)
        if !haskey(st.polygons_dict, a) 
            st.polygons_dict[a] = Set{RationalPolygon{T}}()
        end
        P ∉ st.polygons_dict[a] || continue
        push!(st.polygons_dict[a], P)
        st.total_count += 1
    end
end

is_finished(st :: InMemorySubpolygonStorage{T}) where {T <: Integer} =
st.last_volume <= minimum(keys(st.polygons_dict))


function subpolygons_single_step(
        st :: InMemorySubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    a = maximum(filter(b -> b < st.last_volume, keys(st.polygons_dict)))
    Ps = collect(st.polygons_dict[a])

    logging && @info "[Area: $a]. Number of polygons: $(length(Ps)). Total: $(st.total_count)"

    N = length(Ps)
    num_blocks = 2 * Threads.nthreads()
    b = N ÷ num_blocks

    out_array = Vector{RationalPolygon{T}}[]
    for i = 1 : Threads.nthreads()
        push!(out_array, RationalPolygon{T}[])
    end

    Threads.@threads for k = 1 : num_blocks
        lower_bound = (k-1) * b + 1
        upper_bound = k == num_blocks ? N : k * b
        for i = lower_bound : upper_bound
            P = Ps[i]
            for j = 1 : number_of_vertices(P)
                Q, keeps_genus = remove_vertex(P, j; st.primitive)
                keeps_genus || continue
                number_of_vertices(Q) > 2 || continue
                if st.use_affine_normal_form
                    Q = affine_normal_form(Q)
                else
                    Q = unimodular_normal_form(Q)
                end
                push!(out_array[Threads.threadid()], Q)
            end
        end
    end

    save!(st, vcat(out_array...))

    st.last_volume = a

end

function subpolygons(
        st :: InMemorySubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    while !is_finished(st)
        subpolygons_single_step(st; logging)
    end

    return collect(union!(values(st.polygons_dict)...))

end
