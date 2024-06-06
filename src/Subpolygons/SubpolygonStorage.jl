@doc raw"""
    SubpolygonStorage{T <: Integer}

An abstract supertype of storage options for computing subpolygons. There are
two subtypes `InMemorySubpolygonStorage` and `OnDiskSubpolygonStorage`. The
former keeps all subpolygons in memory, the latter delegeates their storage to
the disk.

"""
abstract type SubpolygonStorage{T <: Integer} end

rationality(st :: SubpolygonStorage{T}) where {T <: Integer} = st.rationality
number_of_interior_lattice_points(st :: SubpolygonStorage{T}) where {T <: Integer} = st.number_of_interior_lattice_points
total_count(st :: SubpolygonStorage{T}) where {T <: Integer} = st.total_count

mutable struct InMemorySubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    rationality :: T
    number_of_interior_lattice_points :: Int
    hashes_dict :: Dict{T,Set{UInt}}
    polygons_dict :: Dict{T,Vector{RationalPolygon{T}}}
    last_volume :: T
    total_count :: Int

    @doc raw"""
        InMemorySubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    Construct an in-memory storage container of subpolygons starting with the
    given list of polygons. When using this constructor, the subpolygons are
    not yet computed. To compute the subpolygons, use [`subpolygons`](@ref).
    
    """
    function InMemorySubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")

        k = rationality(first(starting_polygons))
        n = number_of_interior_lattice_points(first(starting_polygons))
        all(P -> rationality(P) == k, starting_polygons) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, starting_polygons) || error("all polygons must have the same number of interior lattice points")


        hashes_dict = Dict{T, Set{UInt}}()
        polygons_dict = Dict{T, Vector{RationalPolygon{T}}}()
        last_volume = maximum(normalized_area.(starting_polygons)) + 1
        total_count = 0 
        st = new{T}(k, n, hashes_dict, polygons_dict, last_volume, total_count)
        save!(st, starting_polygons)
        return st
    end

end

function save!(st :: InMemorySubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    for P ∈ Ps
        a = normalized_area(P)
        if !haskey(st.hashes_dict, a) 
            st.hashes_dict[a] = Set{UInt}()
            st.polygons_dict[a] = RationalPolygon{T}[]
        end
        h = hash(P)
        h ∉ st.hashes_dict[a] || continue
        push!(st.hashes_dict[a], h)
        push!(st.polygons_dict[a], P)
        st.total_count += 1
    end
end

is_finished(st :: InMemorySubpolygonStorage{T}) where {T <: Integer} =
st.last_volume <= minimum(keys(st.polygons_dict))

function next_polygons(st :: InMemorySubpolygonStorage{T}) where {T <: Integer}
    a = maximum(filter(b -> b < st.last_volume, keys(st.polygons_dict)))
    return (st.polygons_dict[a], a)
end

function mark_volume_completed(st :: InMemorySubpolygonStorage{T}, a :: T) where {T <: Integer}
    st.last_volume = a
end

return_value(st :: InMemorySubpolygonStorage) = vcat(values(st.polygons_dict)...)


mutable struct OnDiskSubpolygonStorage{T<:Integer} <: SubpolygonStorage{T}
    rationality :: T
    number_of_interior_lattice_points :: Int
    hashes_dict :: Dict{T, Set{UInt}}
    directory :: String
    file_prefix :: String
    files :: Dict{T, IOStream}
    last_volume :: T
    total_count :: Int

    @doc raw"""
        OnDiskSubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}, directory :: String; file_prefix :: String = "vol_") where {T <: Integer}

    Construct a new on-disk subpolygon storage container at the given
    directory, starting with the given list of polygons.
    
    """
    function OnDiskSubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}, directory :: String; file_prefix :: String = "vol_") where {T <: Integer}
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")
        isdir(directory) || error("$directory is not a directory")
        isempty(readdir(directory)) || error("$directory is dirty. Please provide an empty directory")

        k = rationality(first(starting_polygons))
        n = number_of_interior_lattice_points(first(starting_polygons))
        all(P -> rationality(P) == k, starting_polygons) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, starting_polygons) || error("all polygons must have the same number of interior lattice points")

        hashes_dict = Dict{T,Set{UInt}}()
        files = Dict{T, IOStream}()
        last_volume = maximum(normalized_area.(starting_polygons)) + 1
        total_count = 0
        st = new{T}(k, n, hashes_dict, directory, file_prefix, files, last_volume, total_count)
        save!(st, starting_polygons)
        return st
    end

    @doc raw"""
        OnDiskSubpolygonStorage{T}(k :: T, n :: T, directory :: String; file_prefix :: String = "vol_") where {T <: Integer}

    This constructor can be used to resume an unfinished computations of
    subpolygons. It takes the rationality `k` and the number of interior
    lattice points `n` as well as the file path where [`subpolygons`](@ref) has
    stored its unfinished results. The resulting object can the be passed to
    [`subpolygons`](@ref) and it will resume the computation where it left off.
    
    """
    function OnDiskSubpolygonStorage{T}(k :: T, n :: T, directory :: String; file_prefix :: String = "vol_") where {T <: Integer}

        isdir(directory) || error("$directory is not a directory")
        
        last_volume_filepath = joinpath(directory, "last_volume.txt")
        isfile(last_volume_filepath) || error("$last_volume_file not found")

        f = open(last_volume_filepath, "r")
        last_volume = parse(Int, readline(f))
        close(f)

        hashes_dict = Dict{T, Set{UInt}}()
        files = Dict{T, IOStream}()
        total_count = 0
        for filename ∈ readdir(directory)
            filename == "last_volume.txt" && continue
            a = parse(Int, replace(filename, file_prefix => "", ".txt" => ""))
            filepath = joinpath(directory, filename)
            total_count += countlines(filepath)
            a < last_volume || continue
            Ps = parse_rational_polygons(k, filepath)
            hashes_dict[a] = Set(hash.(Ps))
            files[a] = open(joinpath(directory, filename), "a")
        end

        return new{T}(k, n, hashes_dict, directory, file_prefix, files, last_volume, total_count)

    end


end


function save!(st :: OnDiskSubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    as = unique(normalized_area.(Ps))
    files = Dict([a => open(joinpath(st.directory, "$(st.file_prefix)$a.txt"), "a") for a ∈ as])

    for P ∈ Ps
        P = unimodular_normal_form(P)
        a = normalized_area(P)
        if !haskey(st.hashes_dict, a) 
            st.hashes_dict[a] = Set{RationalPolygon{T}}()
        end
        h = hash(P)
        h ∉ st.hashes_dict[a] || continue
        push!(st.hashes_dict[a], h)
        V = vertex_matrix(P)
        println(files[a], [Vector(V[:,i]) for i = 1 : number_of_vertices(P)])
        st.total_count += 1
    end

    for a ∈ as
        close(files[a])
    end
end

is_finished(st :: OnDiskSubpolygonStorage{T}) where {T <: Integer} =
isempty(st.hashes_dict)

function next_polygons(st :: OnDiskSubpolygonStorage{T}) where {T <: Integer}
    # Get the polygons with maximal normalized_area
    a = maximum(filter(b -> b < st.last_volume, keys(st.hashes_dict)))

    filepath = joinpath(st.directory, "$(st.file_prefix)$a.txt")
    Ps = parse_rational_polygons(st.rationality, filepath)

    return (Ps, a)

end

function mark_volume_completed(st :: OnDiskSubpolygonStorage{T}, a :: T) where {T <: Integer}
    delete!(st.hashes_dict, a)
    st.last_volume = a

    open(joinpath(st.directory, "last_volume.txt"), "w") do f
        println(f, st.last_volume)
    end
end

return_value(st :: OnDiskSubpolygonStorage) = nothing

subpolygon_storage(starting_polygons :: Vector{<:RationalPolygon{T}}, directory :: Union{Missing,String} = missing) where {T <: Integer} =
ismissing(directory) ? InMemorySubpolygonStorage{T}(starting_polygons) : OnDiskSubpolygonStorage{T}(starting_polygons, directory)
