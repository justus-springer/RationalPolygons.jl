@doc raw"""
    parse_rational_polygons_file(k :: T, filepath :: String) where {T <: Integer}

Parse a file containing the vertices of a `k`-rational polygon. The file
must have one polygon per line and the vertices must be given as a list of
lists of integers, i.e. as in the following example:

[[2, 0], [1, 3], [-1, 0], [-3, -4]]
[[1, 0], [2, 6], [-4, -9]]
[[1, 0], [3, 5], [0, 1], [-5, -8]]
....

"""
function parse_rational_polygons_file(k :: T, filepath :: String) where {T <: Integer}
    f = open(filepath, "r")
    Ps = RationalPolygon{T}[]
    for str ∈ readlines(f)
        substrings = map(s -> replace(s, "[" => "", "]" => ""), split(str, ", ["))
        vs = [parse.(Int,split(s, ",")) for s ∈ substrings]
        P = RationalPolygon([LatticePoint{T}(v[1],v[2]) for v ∈ vs], k)
        push!(Ps, P)
    end
    close(f)
    return Ps
end


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

    function InMemorySubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")

        k = rationality(first(starting_polygons))
        n = number_of_interior_lattice_points(first(starting_polygons))
        all(P -> rationality(P) == k, starting_polygons) || error("all polygons must have the same rationality")
        all(P -> number_of_interior_lattice_points(P) == n, starting_polygons) || error("all polygons must have the same number of interior lattice points")


        hashes_dict = Dict{T, Set{UInt}}()
        polygons_dict = Dict{T, Vector{RationalPolygon{T}}}()
        last_volume = maximum(area.(starting_polygons)) + 1
        total_count = 0 
        st = new{T}(k, n, hashes_dict, polygons_dict, last_volume, total_count)
        save!(st, starting_polygons)
        return st
    end

end

function save!(st :: InMemorySubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    for P ∈ Ps
        a = area(P)
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

function next_polygons!(st :: InMemorySubpolygonStorage{T}) where {T <: Integer}
    a = maximum(filter(b -> b < st.last_volume, keys(st.polygons_dict)))
    return (st.polygons_dict[a], a)
end

function mark_volume_completed(st :: InMemorySubpolygonStorage{T}, a :: T) where {T <: Integer}
    st.last_volume = a
end

return_value(st :: InMemorySubpolygonStorage) = vcat(values(st.polygons_dict)...)


@doc raw"""
    mutable struct OnDiskSubpolygonStorage{T<:Integer}   

A struct used by the function `subpolygons` to hold intermediate results.

"""
mutable struct OnDiskSubpolygonStorage{T<:Integer} <: SubpolygonStorage{T}
    rationality :: T
    number_of_interior_lattice_points :: Int
    hashes_dict :: Dict{T, Set{UInt}}
    directory :: String
    file_prefix :: String
    files :: Dict{T, IOStream}
    last_volume :: T
    total_count :: Int

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
        last_volume = maximum(area.(starting_polygons)) + 1
        total_count = 0
        st = new{T}(k, n, hashes_dict, directory, file_prefix, files, last_volume, total_count)
        save!(st, starting_polygons)
        return st
    end

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
            Ps = parse_rational_polygons_file(k, filepath)
            hashes_dict[a] = Set(hash.(Ps))
            files[a] = open(joinpath(directory, filename), "a")
        end

        return new{T}(k, n, hashes_dict, directory, file_prefix, files, last_volume, total_count)

    end


end


function save!(st :: OnDiskSubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    as = unique(area.(Ps))
    files = Dict([a => open(joinpath(st.directory, "$(st.file_prefix)$a.txt"), "a") for a ∈ as])

    for P ∈ Ps
        P = normal_form(P)
        a = area(P)
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

function next_polygons(st :: OnDiskSubpolygonStorage{T}) where {T <: Integer}
    # Get the polygons with maximal area
    a = maximum(filter(b -> b < st.last_volume, keys(st.hashes_dict)))

    filepath = joinpath(st.directory, "$(st.file_prefix)$a.txt")
    Ps = parse_rational_polygons_file(st.rationality, filepath)

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
