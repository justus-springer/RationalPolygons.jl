
@doc raw"""
    parse_rational_polygons(k :: T, files :: AbstractVector{String}) where {T <: Integer}
    parse_rational_polygons(k :: T, file :: String) where {T <: Integer}

Parse a list of files containing the vertices of a `k`-rational polygon. The
files must have one polygon per line and the vertices must be given as a list of
lists of integers, i.e. as in the following example:

[[2, 0], [1, 3], [-1, 0], [-3, -4]]
[[1, 0], [2, 6], [-4, -9]]
[[1, 0], [3, 5], [0, 1], [-5, -8]]
....

"""
function parse_rational_polygons(k :: T, files :: AbstractVector{String}) where {T <: Integer}
    Ps = RationalPolygon{T}[]
    for fp ∈ files
        f = open(fp, "r")
        for str ∈ readlines(f)
            substrings = map(s -> replace(s, "[" => "", "]" => ""), split(str, ", ["))
            vs = [parse.(Int,split(s, ",")) for s ∈ substrings]
            P = RationalPolygon([LatticePoint{T}(v[1],v[2]) for v ∈ vs], k)
            push!(Ps, P)
        end
        close(f)
    end
    return Ps
end

parse_rational_polygons(k :: T, file :: String) where {T <: Integer} =
parse_rational_polygons(k, [file])


@doc raw"""
    write_rational_polygons(Ps :: Vector{<:RationalPolygon{T}}, filepath :: String) where {T <: Integer}

Write a list of polygons to a text file. Each polygon will be written as one
line containing its vertices, as in the following example:

[[2, 0], [1, 3], [-1, 0], [-3, -4]]
[[1, 0], [2, 6], [-4, -9]]
[[1, 0], [3, 5], [0, 1], [-5, -8]]
....

"""
function write_rational_polygons(Ps :: Vector{<:RationalPolygon{T}}, filepath :: String) where {T <: Integer}
   f = open(filepath, "w")
   for P ∈ Ps  
       V = vertex_matrix(P)
       println(f, [Vector(V[:,i]) for i = 1 : number_of_vertices(P)])
   end
   close(f)
end
