
abstract type SubpolygonStorage{T <: Integer} end

mutable struct InMemorySubpolygonStorage{T <: Integer} <: SubpolygonStorage{T}
    hashes_dict :: Dict{T,Set{UInt}}
    polygons_dict :: Dict{T,Vector{RationalPolygon{T}}}
    last_volume :: T

    function InMemorySubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")

        a = maximum(area.(starting_polygons))
        st = new{T}(Dict{T, Set{UInt64}}(), Dict{T, Vector{RationalPolygon{T}}}(), a + 1)
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
    end
end

function next_polygons!(st :: InMemorySubpolygonStorage{T}) where {T <: Integer}
    a = maximum(filter(b -> b < st.last_volume, keys(st.polygons_dict)))
    st.last_volume = a
    return (st.polygons_dict[a], a)
end


@doc raw"""
    mutable struct OnDiskSubpolygonStorage{T<:Integer}   

A struct used by the function `subpolygons` to hold intermediate results.

"""
mutable struct OnDiskSubpolygonStorage{T<:Integer} <: SubpolygonStorage{T}
    rationality :: T
    hashes_dict :: Dict{T, Set{UInt}}
    directory :: String
    file_prefix :: String
    files :: Dict{T, IOStream}
    last_volume :: T

    function OnDiskSubpolygonStorage{T}(starting_polygons :: Vector{<:RationalPolygon{T}}, directory :: String; file_prefix :: String = "vol_") where {T <: Integer}
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")
        isdir(directory) || error("$directory is not a directory")

        a = maximum(area.(starting_polygons))
        k = rationality(first(starting_polygons))
        st = new{T}(k, Dict{T,Set{UInt}}(), directory, file_prefix, Dict{T, IOStream}(), a + 1)
        save!(st, starting_polygons)
        return st
    end

end

function save!(st :: OnDiskSubpolygonStorage{T}, Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    for P ∈ Ps
        P = normal_form(P)
        a = area(P)
        if !haskey(st.hashes_dict, a) 
            st.hashes_dict[a] = Set{RationalPolygon{T}}()
            st.files[a] = open(joinpath(st.directory, "$(st.file_prefix)$a.txt"), "a")
        end
        h = hash(P)
        h ∉ st.hashes_dict[a] || continue
        push!(st.hashes_dict[a], h)
        V = vertex_matrix(P)
        println(st.files[a], [Vector(V[:,i]) for i = 1 : number_of_vertices(P)])
    end
end

function next_polygons!(st :: SubpolygonStorage{T}) where {T <: Integer}
    # Get the polygons with maximal area
    a = maximum(filter(b -> b < st.last_volume, keys(st.hashes_dict)))
    close(st.files[a])

    f = open(joinpath(st.directory, "$(st.file_prefix)$a.txt"), "r")
    Ps = RationalPolygon{T}[]
    for str ∈ readlines(f)
        substrings = map(s -> replace(s, "[" => "", "]" => ""), split(str, ", ["))
        vs = [parse.(Int,split(s, ",")) for s ∈ substrings]
        P = RationalPolygon([LatticePoint{T}(v[1],v[2]) for v ∈ vs], st.rationality)
        push!(Ps, P)
    end
    close(f)

    delete!(st.hashes_dict, a)
    delete!(st.files, a)
    st.last_volume = a
    return (Ps, a)

end

subpolygon_storage(starting_polygons :: Vector{<:RationalPolygon{T}}, directory :: Union{Missing,String} = missing) where {T <: Integer} =
ismissing(directory) ? InMemorySubpolygonStorage{T}(starting_polygons) : OnDiskSubpolygonStorage{T}(starting_polygons, directory)
