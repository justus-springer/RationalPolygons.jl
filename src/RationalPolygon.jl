
@attributes mutable struct RationalPolygon{T<:Integer}
    rationality :: T
    vertex_matrix :: Matrix{T}

    RationalPolygon(vertex_matrix :: Matrix{T}, rationality :: T) where {T <: Integer} =
    new{T}(rationality, vertex_matrix)

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

rationality(P :: RationalPolygon) = P.rationality

vertex_matrix(P :: RationalPolygon) = P.vertex_matrix

function Base.hash(P :: RationalPolygon, h :: UInt64)
    h = hash(rationality(P), h)
    h = hash(vertex_matrix(P), h)
    return hash(RationalPolygon, h)
end

Base.:(==)(P1 :: RationalPolygon, P2 :: RationalPolygon) =
rationality(P1) == rationality(P2) && vertex_matrix(P1) == vertex_matrix(P2)

Base.show(io :: IO, P :: RationalPolygon) =
Base.print(io, "Rational polygon of rationality $(rationality(P)) with $(number_of_vertices(P)) vertices.")

@attr number_of_vertices(P :: RationalPolygon) = size(vertex_matrix(P), 2)

function vertex(P :: RationalPolygon, i :: Int)
    V, n, k = vertex_matrix(P), number_of_vertices(P), rationality(P)
    i = mod(i, 1:n)
    return (V[1,i] // k, V[2,i] // k)
end

Base.getindex(P :: RationalPolygon, i :: Int) = vertex(P, i)

Base.iterate(P :: RationalPolygon) = (P[1], 2)
function Base.iterate(P :: RationalPolygon, i :: Int) 
    if i > number_of_vertices(P) 
        return nothing
    else
        return (P[i], i+1)
    end
end

Base.IteratorSize(:: Type{<:RationalPolygon}) = Base.HasLength()
Base.length(P :: RationalPolygon) = number_of_vertices(P)

Base.IteratorEltype(:: Type{<:RationalPolygon}) = Base.HasEltype()
Base.eltype(:: Type{<:RationalPolygon{T}}) where {T <: Integer} =
RationalPoint{T}

vertices(P :: RationalPolygon) = collect(P)

@attr function affine_halfplanes(P :: RationalPolygon)
    n = number_of_vertices(P)
    return [AffineHalfplaneByLine(LineThroughPoints(P[i], P[mod(i+1,1:n)])) for i = 1 : n]
end

Base.in(x :: Point{T}, P :: RationalPolygon) where {T <: Integer} =
all(H -> x ∈ H, affine_halfplanes(P))

@attr function is_primitive(P :: RationalPolygon)
    V = vertex_matrix(P)
    return all(i -> gcd(V[1,i], V[2,i]) == 1, 1 : number_of_vertices(P))
end

@doc raw"""
    is_maximal(P :: RationalPolygon)

Check whether a rational polygon is maximal among all polygons sharing
the same rationality and number of interior lattice points.

"""
@attr function is_maximal(P :: RationalPolygon)
    k = rationality(P)
    ps = k_rational_points(k,P)
    Q = intersect_halfplanes(affine_halfplanes(P) .- 1 // k)
    for p ∈ boundary_k_rational_points(k, Q)
        new_P = convex_hull([ps ; p], k)
        number_of_interior_lattice_points(new_P) <= number_of_interior_lattice_points(P) && return false
    end
    return true
end
