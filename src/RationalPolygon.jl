
struct RationalPolygon{T<:Integer,N,M}
    rationality :: T
    vertex_matrix :: SMatrix{2, N, T, M}
    is_normal_form :: Bool

    RationalPolygon(vertex_matrix :: SMatrix{2,N,T,M},
                    rationality :: T; 
                    is_normal_form :: Bool = false) where {N, M, T <: Integer} =
    new{T,N,M}(rationality, vertex_matrix, is_normal_form)

    RationalPolygon(scaled_points :: Vector{LatticePoint{T}}, 
                    rationality :: T; 
                    is_normal_form :: Bool = false) where {T <: Integer} =
    RationalPolygon(hcat(scaled_points...), rationality)

    function RationalPolygon(points :: Vector{RationalPoint{T}};
            is_normal_form :: Bool = false) where {T <: Integer}
        k = lcm(rationality.(points))
        scaled_points = numerator.(k .* points)
        return RationalPolygon(scaled_points, k; is_normal_form)
    end

    function RationalPolygon(points :: Vector{RationalPoint{T}}, 
            rationality :: T;
            is_normal_form :: Bool = false) where {T <: Integer}
        scaled_points = numerator.(rationality .* points)
        return RationalPolygon(scaled_points, rationality; is_normal_form)
    end

end

empty_polygon(rationality :: T) where {T <: Integer} =
RationalPolygon(SMatrix{2,0,T,0}(), rationality)

empty_polygon(::Type{T}) where {T <: Integer} =
empty_polygon(one(T))

convex_hull(points :: Vector{LatticePoint{T}}, k :: T) where {T <: Integer} =
RationalPolygon(graham_scan(points), k)

convex_hull(points :: Vector{RationalPoint{T}}) where {T <: Integer} =
RationalPolygon(graham_scan(points))

convex_hull(points :: Vector{RationalPoint{T}}, k :: T) where {T <: Integer} =
RationalPolygon(graham_scan(points), k)

number_of_vertices(P :: RationalPolygon{T,N,M}) where {N,M,T <: Integer} = N

rationality(P :: RationalPolygon{T,N}) where {N,T <: Integer} = P.rationality

vertex_matrix(P :: RationalPolygon{T,N}) where {N,T <: Integer} = P.vertex_matrix

is_normal_form(P :: RationalPolygon) = P.is_normal_form

clockwise(P :: RationalPolygon{T,N}) where {N, T <: Integer} = P.clockwise

Base.:(==)(P1 :: RationalPolygon{T,N}, P2 :: RationalPolygon{T,N}) where {N,T <: Integer} =
rationality(P1) == rationality(P2) && vertex_matrix(P1) == vertex_matrix(P2)

Base.show(io :: IO, P :: RationalPolygon{T,N}) where {N,T <: Integer} =
Base.print(io, "Rational polygon of rationality $(rationality(P)) with $(number_of_vertices(P)) vertices.")

function lattice_vertex(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer}
    V = vertex_matrix(P)
    i = mod(i, 1:N)
    return V[:,i]
end

vertex(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} = lattice_vertex(P, i) .// rationality(P)

function Base.hash(P :: RationalPolygon{T,N}, h :: UInt64) where {N,T <: Integer}
    h = hash(P.rationality, h)
    for i = 1 : N
        h = hash(lattice_vertex(P,i), h)
    end
    return hash(RationalPolygon, h)
end

Base.getindex(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} = vertex(P, i)

Base.iterate(P :: RationalPolygon{T,N}) where {N,T <: Integer} = (P[1], 2)
Base.iterate(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} =
i > N ? nothing : (P[i], i+1)

Base.IteratorSize(:: Type{<:RationalPolygon{T,N}}) where {N,T <: Integer} = Base.HasLength()
Base.length(P :: RationalPolygon{T,N}) where {N,T <: Integer} = N

Base.IteratorEltype(:: Type{<:RationalPolygon{T,N}}) where {N,T <: Integer} = Base.HasEltype()
Base.eltype(:: Type{<:RationalPolygon{T,N}}) where {N,T <: Integer} =
RationalPoint{T}

vertices(P :: RationalPolygon{T,N}) where {N,T <: Integer} = collect(P)

affine_halfplanes(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
[affine_halfplane(line_through_points(P[i], P[mod(i+1,1:N)])) for i = 1 : N]

Base.in(x :: Point{T}, P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(H -> x ∈ H, affine_halfplanes(P))

is_primitive(P :: RationalPolygon{T,N}) where {N,T <: Integer} =
all(i -> is_primitive(lattice_vertex(P,i)), 1 : N)

function area(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    V = vertex_matrix(P)
    return sum([abs(det(V[:,i+1] - V[:,1], V[:,i] - V[:,1])) for i = 2 : N -1])
end

@doc raw"""
    is_maximal(P :: RationalPolygon)

Check whether a rational polygon is maximal among all polygons sharing
the same rationality and number of interior lattice points.

"""
function is_maximal(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    k = rationality(P)
    ps = k_rational_points(k,P)
    Q = intersect_halfplanes(affine_halfplanes(P) .- 1 // k)
    for p ∈ boundary_k_rational_points(k, Q)
        new_P = convex_hull([ps ; [p]], k)
        number_of_interior_lattice_points(new_P) <= number_of_interior_lattice_points(P) && return false
    end
    return true
end
