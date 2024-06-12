
mutable struct OnDiskSubpolygonStorage{T<:Integer} <: SubpolygonStorage{T}
    rationality :: T
    number_of_interior_lattice_points :: Int
    primitive :: Bool
    normal_form :: Symbol
    directory :: String
    active_volumes :: Vector{T}
    total_count :: Int

    @doc raw"""
        OnDiskSubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}, directory :: String)

    Construct a new on-disk subpolygon storage container at the given
    directory, starting with the given list of polygons.
    
    """
    function OnDiskSubpolygonStorage{T}(
            starting_polygons :: Vector{<:RationalPolygon{T}}, 
            directory :: String;
            primitive :: Bool = false,
            normal_form :: Symbol = :unimodular) where {T <: Integer}

        normal_form ∈ [:unimodular, :affine] || error("`normal_form` must be either `:unimodular` or `:affine`")
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")
        isdir(directory) || error("$directory is not a directory")
        isempty(readdir(directory)) || error("$directory is dirty. Please provide an empty directory")

        k = rationality(first(starting_polygons))
        n = number_of_interior_lattice_points(first(starting_polygons))
        all(P -> rationality(P) == k, starting_polygons) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, starting_polygons) || error("all polygons must have the same number of interior lattice points")

        st = new{T}(k, n, primitive, normal_form, directory, T[], 0)
        save!(st, starting_polygons)
        return st
    end

    #=
    @doc raw"""
        OnDiskSubpolygonStorage{T}(k :: T, n :: T, directory :: String; file_prefix :: String = "vol_") where {T <: Integer}

    This constructor can be used to resume an unfinished computations of
    subpolygons. It takes the rationality `k` and the number of interior
    lattice points `n` as well as the file path where [`subpolygons`](@ref) has
    stored its unfinished results. The resulting object can the be passed to
    [`subpolygons`](@ref) and it will resume the computation where it left off.
    
    """
    function OnDiskSubpolygonStorage{T}(
            k :: T, 
            n :: T, 
            directory :: String; 
            primitive :: Bool = false,
            normal_form :: Symbol = :unimodular) where {T <: Integer}

        isdir(directory) || error("$directory is not a directory")
        
        last_volume_filepath = joinpath(directory, "last_volume.txt")
        isfile(last_volume_filepath) || error("$last_volume_file not found")
        f = open(last_volume_filepath, "r")
        last_volume = parse(Int, readline(f))
        close(f)

        total_count_filepath = joinpath(directory, "total_count.txt")
        isfile(total_count_filepath) || error("$total_count_filepath not found")
        f = open(total_count_filepath, "r")
        total_count = parse(Int, readline(f))
        close(f)

        active_volumes = T[]
        for filename ∈ readdir(directory)
            startswith(filename, "vol_") || continue
            a = parse(T, replace(filename, file_prefix => "", ".txt" => ""))
            filepath = joinpath(directory, filename)
            total_count += countlines(filepath)
            a < last_volume || continue
            Ps = parse_rational_polygons(k, filepath)
            hashes_dict[a] = Set(hash.(Ps))
            files[a] = open(joinpath(directory, filename), "a")
        end

        return new{T}(k, n, hashes_dict, directory, file_prefix, files, last_volume, total_count)

    end
    =#

end

function save!(st :: OnDiskSubpolygonStorage{T}, dict :: Dict{T, Set{RationalPolygon{T}}}) where {T <: Integer}
    for (a,Ps) ∈ dict
        filepath = joinpath(st.directory, "vol_$a.txt")
        if isfile(filepath)
            Qs = Set{RationalPolygon{T}}(parse_rational_polygons(rationality(st), filepath))
        else
            Qs = Set{RationalPolygon{T}}()
            push!(st.active_volumes, a)
        end

        f = open(filepath, "a")
        for P ∈ Ps
            P ∉ Qs || continue
            V = vertex_matrix(P)
            println(f, [Vector(V[:,i]) for i = 1 : number_of_vertices(P)])
            st.total_count += 1
        end
        close(f)

    end
end

function save!(st :: OnDiskSubpolygonStorage{T}, Ps :: Vector{RationalPolygon{T}}) where {T <: Integer}
    dict = Dict{T,Set{RationalPolygon{T}}}()
    for P ∈ Ps
        a = normalized_area(P)
        if !haskey(dict, a)
            dict[a] = Set{RationalPolygon{T}}()
        end
        push!(dict[a], P)
    end
    save!(st, dict)
end

is_finished(st :: OnDiskSubpolygonStorage{T}) where {T <: Integer} =
isempty(st.active_volumes)


function subpolygons_single_step(
        st :: OnDiskSubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    a = maximum(st.active_volumes)
    filepath = joinpath(st.directory, "vol_$a.txt")
    Ps = parse_rational_polygons(st.rationality, filepath)

    logging && @info "[Area: $a]. Number of polygons: $(length(Ps)). Total: $(total_count(st))"

    N = length(Ps)
    num_blocks = 2 * Threads.nthreads()
    b = N ÷ num_blocks

    dicts = Dict{T,Set{RationalPolygon{T}}}[]
    for i = 1 : Threads.nthreads()
        push!(dicts, Dict{T,Set{RationalPolygon{T}}}())
    end

    Threads.@threads for k = 1 : num_blocks
        id = Threads.threadid()
        lower_bound = (k-1) * b + 1
        upper_bound = k == num_blocks ? N : k * b
        for i = lower_bound : upper_bound
            P = Ps[i]
            for j = 1 : number_of_vertices(P)
                Q, keeps_genus = remove_vertex(P, j; st.primitive)
                keeps_genus || continue
                number_of_vertices(Q) > 2 || continue
                if st.normal_form == :unimodular
                    Q = unimodular_normal_form(Q)
                else
                    Q = affine_normal_form(Q)
                end

                area = normalized_area(Q)
                if !haskey(dicts[id], area) 
                    dicts[id][area] = Set{RationalPolygon{T}}()
                end
                Q ∉ dicts[id][area] || continue
                push!(dicts[id][area], Q)
            end
        end
    end
    total_dict = mergewith(union!, dicts...)
    save!(st, total_dict)

    filter!(b -> b < a, st.active_volumes)
    open(joinpath(st.directory, "last_volume.txt"), "w") do f
        println(f, a)
    end
    open(joinpath(st.directory, "total_count.txt"), "w") do f
        println(f, st.total_count)
    end

end


return_value(st :: OnDiskSubpolygonStorage) = nothing
