
@attributes mutable struct RationalPolygon{T<:Integer}
    rationality :: T
    vertex_matrix :: Matrix{T}
    vertex_offset :: Int
    clockwise :: Bool

    RationalPolygon(vertex_matrix :: Matrix{T}, rationality :: T; vertex_offset :: Int = 1, clockwise :: Bool = false) where {T <: Integer} =
    new{T}(rationality, vertex_matrix, vertex_offset, clockwise)

    function RationalPolygon(scaled_points :: Vector{LatticePoint{T}}, rationality :: T) where {T <: Integer}
        n = length(scaled_points)
        V = Matrix{T}(undef, 2, n)
        for i = 1 : n
            V[1,i] = scaled_points[i][1]
            V[2,i] = scaled_points[i][2]
        end
        return RationalPolygon(V, rationality)
    end

    function RationalPolygon(points :: Vector{RationalPoint{T}}) where {T <: Integer}
        k = lcm(rationality.(points))
        scaled_points = numerator.(k .* points)
        return RationalPolygon(scaled_points, k)
    end

    function RationalPolygon(points :: Vector{RationalPoint{T}}, rationality :: T) where {T <: Integer}
        scaled_points = numerator.(rationality .* points)
        return RationalPolygon(scaled_points, rationality)
    end

end

empty_polygon(rationality :: T) where {T <: Integer} =
RationalPolygon(LatticePoint{T}[], rationality)

empty_polygon(::Type{T}) where {T <: Integer} =
RationalPolygon(LatticePoint{T}[], 1)

convex_hull(points :: Vector{LatticePoint{T}}, k :: T) where {T <: Integer} =
RationalPolygon(graham_scan(points), k)

convex_hull(points :: Vector{RationalPoint{T}}) where {T <: Integer} =
RationalPolygon(graham_scan(points))

convex_hull(points :: Vector{RationalPoint{T}}, k :: T) where {T <: Integer} =
RationalPolygon(graham_scan(points), k)

rationality(P :: RationalPolygon{T}) where {T <: Integer} = P.rationality

vertex_matrix(P :: RationalPolygon{T}) where {T <: Integer} = P.vertex_matrix

vertex_offset(P :: RationalPolygon{T}) where {T <: Integer} = P.vertex_offset

clockwise(P :: RationalPolygon{T}) where {T <: Integer} = P.clockwise

Base.:(==)(P1 :: RationalPolygon{T}, P2 :: RationalPolygon{T}) where {T <: Integer} =
rationality(P1) == rationality(P2) && vertex_matrix(P1) == vertex_matrix(P2)

Base.show(io :: IO, P :: RationalPolygon{T}) where {T <: Integer} =
Base.print(io, "Rational polygon of rationality $(rationality(P)) with $(number_of_vertices(P)) vertices.")

number_of_vertices(P :: RationalPolygon{T}) where {T <: Integer} = size(vertex_matrix(P), 2)

function lattice_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
    V, n = vertex_matrix(P), number_of_vertices(P)
    o = vertex_offset(P)
    s = clockwise(P) ? -1 : 1
    i = mod(o + s * (i-1), 1:n)
    return Tuple{T,T}((V[1,i], V[2,i]))
end

vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} = lattice_vertex(P, i) .// rationality(P)

function Base.hash(P :: RationalPolygon{T}, h :: UInt64) where {T <: Integer}
    h = hash(P.rationality, h)
    for i = 1 : number_of_vertices(P)
        h = hash(lattice_vertex(P,i), h)
    end
    return hash(RationalPolygon, h)
end

Base.getindex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer} = vertex(P, i)

Base.iterate(P :: RationalPolygon{T}) where {T <: Integer} = (P[1], 2)
function Base.iterate(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
    if i > number_of_vertices(P) 
        return nothing
    else
        return (P[i], i+1)
    end
end

Base.IteratorSize(:: Type{<:RationalPolygon{T}}) where {T <: Integer} = Base.HasLength()
Base.length(P :: RationalPolygon{T}) where {T <: Integer} = number_of_vertices(P)

Base.IteratorEltype(:: Type{<:RationalPolygon{T}}) where {T <: Integer} = Base.HasEltype()
Base.eltype(:: Type{<:RationalPolygon{T}}) where {T <: Integer} =
RationalPoint{T}

vertices(P :: RationalPolygon{T}) where {T <: Integer} = collect(P)

@attr function affine_halfplanes(P :: RationalPolygon{T}) where {T <: Integer}
    n = number_of_vertices(P)
    return [AffineHalfplaneByLine(LineThroughPoints(P[i], P[mod(i+1,1:n)])) for i = 1 : n]
end

Base.in(x :: Point{T}, P :: RationalPolygon{T}) where {T <: Integer} =
all(H -> x ∈ H, affine_halfplanes(P))

@attr function is_primitive(P :: RationalPolygon{T}) where {T <: Integer}
    V = vertex_matrix(P)
    return all(i -> gcd(V[1,i], V[2,i]) == 1, 1 : number_of_vertices(P))
end

@attr function area(P :: RationalPolygon{T}) where {T <: Integer}
    V, n = vertex_matrix(P), number_of_vertices(P)
    return sum([abs(det((V[1,i+1] - V[1,1], V[2,i+1] - V[2,1]), (V[1,i] - V[1,1], V[2,i] - V[2,1]))) for i = 2 : n -1])
end

@doc raw"""
    is_maximal(P :: RationalPolygon)

Check whether a rational polygon is maximal among all polygons sharing
the same rationality and number of interior lattice points.

"""
@attr function is_maximal(P :: RationalPolygon{T}) where {T <: Integer}
    k = rationality(P)
    ps = k_rational_points(k,P)
    Q = intersect_halfplanes(affine_halfplanes(P) .- 1 // k)
    for p ∈ boundary_k_rational_points(k, Q)
        new_P = convex_hull([ps ; p], k)
        number_of_interior_lattice_points(new_P) <= number_of_interior_lattice_points(P) && return false
    end
    return true
end
