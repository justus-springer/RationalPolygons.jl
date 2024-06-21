
@doc raw"""
    RationalPolygon{T<:Integer,N,M}

The type of rational polygons in two-dimensional space with `N` vertices. The
parameter `M` is always equal to `2*N`.

"""
struct RationalPolygon{T<:Integer,N,M}

    rationality :: T

    number_of_vertices :: Int

    vertex_matrix :: SMatrix{2, N, T, M}

    is_unimodular_normal_form :: Bool

    is_affine_normal_form :: Bool


    @doc raw"""
        RationalPolygon(vertex_matrix :: SMatrix{2,N,T,M}, rationality :: T) where {N, M, T <: Integer}
        RationalPolygon(scaled_vertices :: Vector{LatticePoint{T}}, rationality :: T)
        RationalPolygon(vertices :: Vector{RationalPoint{T}})
        RationalPolygon(vertices :: Vector{RationalPoint{T}}, rationality :: T)

    Construct a rational polygon from a given set of vertices and possibly
    prescribed rationality. If `vertex_matrix` or `scaled_vertices` is
    provided, these are understood to be the integral vertices of the scaled
    polygon `rationality(P) * P`.

    When using this constructor, no consistency checks are done on the input.
    In particular, the user must be sure that the given points are truly
    vertices of the polygon and *that they are in the correct order* (both
    clockwise and counterclockwise is okay). If this is not known ahead of
    contruction, `convex_hull` should be used instead of this constructor.

    All constructors accept the optional arguments `is_unimodular_normal_form`
    and `is_affine_normal_form`, which are set to `false` by default. If they
    are set to true, this means the user is certain that the given polygon is
    already in the respective normal form. This information will be used to
    prevent addional computations of the normal form and thus speed up
    equivalence checking.

    """
    RationalPolygon(vertex_matrix :: SMatrix{2,N,T,M},
                    rationality :: T; 
                    is_unimodular_normal_form :: Bool = false,
                    is_affine_normal_form :: Bool = false) where {N, M, T <: Integer} =
    new{T,N,M}(rationality, N, vertex_matrix, is_unimodular_normal_form, is_affine_normal_form)

    RationalPolygon(scaled_vertices :: Vector{LatticePoint{T}}, 
                    rationality :: T; 
                    is_unimodular_normal_form :: Bool = false,
                    is_affine_normal_form :: Bool = false) where {T <: Integer} =
    RationalPolygon(hcat(scaled_vertices...), rationality; is_unimodular_normal_form, is_affine_normal_form)

    function RationalPolygon(vertices :: Vector{RationalPoint{T}};
            is_unimodular_normal_form :: Bool = false,
            is_affine_normal_form :: Bool = false) where {T <: Integer}
        k = lcm(rationality.(vertices))
        scaled_vertices = numerator.(k .* vertices)
        return RationalPolygon(scaled_vertices, k; is_unimodular_normal_form, is_affine_normal_form)
    end

    function RationalPolygon(vertices :: Vector{RationalPoint{T}}, 
            rationality :: T;
            is_unimodular_normal_form :: Bool = false,
            is_affine_normal_form :: Bool = false) where {T <: Integer}
        scaled_vertices = numerator.(rationality .* vertices)
        return RationalPolygon(scaled_vertices, rationality; is_unimodular_normal_form, is_affine_normal_form)
    end

end


@doc raw"""
    empty_polygon(rationality :: T) where {T <: Integer}

Return the empty polygon of given rationality.

"""
empty_polygon(rationality :: T) where {T <: Integer} =
RationalPolygon(SMatrix{2,0,T,0}(), rationality)


@doc raw"""
empty_polygon(::Type{T}) where {T <: Integer}

Return the empty polygon of integer type `T`. The rationality is understood to
be one, i.e. it is a lattice polygon.

"""
empty_polygon(::Type{T}) where {T <: Integer} =
empty_polygon(one(T))


@doc raw"""
    convex_hull(points :: Vector{LatticePoint{T}}, k :: T) where {T <: Integer}

Return the `k`-rational polygon given by the convex hull of `p // k`, where `p
∈ points`.

"""
convex_hull(points :: Vector{LatticePoint{T}}, k :: T) where {T <: Integer} =
RationalPolygon(graham_scan(points), k)


@doc raw"""
    convex_hull(points :: Vector{RationalPoint{T}}) where {T <: Integer}

Return the convex hull of a given set of rational points. The rationality will
be inferred from the input.

"""
convex_hull(points :: Vector{RationalPoint{T}}) where {T <: Integer} =
RationalPolygon(graham_scan(points))


@doc raw"""
    convex_hull(points :: Vector{RationalPoint{T}}, k :: T) where {T <: Integer}

Return the convex hull of a given set of rational points, viewed as a
`k`-rational polygon.

"""
convex_hull(points :: Vector{RationalPoint{T}}, k :: T) where {T <: Integer} =
RationalPolygon(graham_scan(points), k)


@doc raw"""
    number_of_vertices(P :: RationalPolygon)

Return the number of vertices of a polygon.

"""
number_of_vertices(P :: RationalPolygon{T,N,M}) where {N,M,T <: Integer} = P.number_of_vertices


@doc raw"""
    rationality(P :: RationalPolygon)

Return the rationality of `P`. Note that this does not need to be the smallest
positive integer `k` such that `k*P` is a lattice polygon: The standard lattice
triangle may also be viewed as a half-integral polygon, in which case the
rationality would be `2`, even though all vertices are integral.

"""
rationality(P :: RationalPolygon{T,N}) where {N,T <: Integer} = P.rationality


@doc raw"""
    vertex_matrix(P :: RationalPolygon)

The vertex matrix of `P` is the 2xN integral matrix containing the vertices of
`rationality(P) * P` as its columns.

"""
vertex_matrix(P :: RationalPolygon{T,N}) where {N,T <: Integer} = P.vertex_matrix


@doc raw"""
    is_unimodular_normal_form(P :: RationalPolygon)

Return whether it is already known that `P` is in unimodular normal form.

"""
is_unimodular_normal_form(P :: RationalPolygon) = P.is_unimodular_normal_form


@doc raw"""
    is_affine_normal_form(P :: RationalPolygon)

Return whether it is already known that `P` is in affine normal form.

"""
is_affine_normal_form(P :: RationalPolygon) = P.is_affine_normal_form

Base.:(==)(P1 :: RationalPolygon{T,N}, P2 :: RationalPolygon{T,N}) where {N,T <: Integer} =
rationality(P1) == rationality(P2) && vertex_matrix(P1) == vertex_matrix(P2)

Base.show(io :: IO, P :: RationalPolygon{T,N}) where {N,T <: Integer} =
Base.print(io, "Rational polygon of rationality $(rationality(P)) with $(number_of_vertices(P)) vertices.")


@doc raw"""
    scaled_vertex(P :: RationalPolygon, i :: Int)

Return the `i`-th integral vertex of the polygon `rationality(P) * P`. The
index `i` is regarded as a cyclic index, e.g. the `N+1`-th vertex is equal to
the first vertex.

"""
function scaled_vertex(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer}
    V = vertex_matrix(P)
    i = mod(i, 1:N)
    return V[:,i]
end


@doc raw"""
    vertex(P :: RationalPolygon, i :: Int)

Return the `i`-th vertex of the polygon `P`. The index `i` is regarded as a
cyclic index, e.g. the `N+1`-th vertex is equal to the first vertex.

"""
vertex(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer} = scaled_vertex(P, i) .// rationality(P)

function Base.hash(P :: RationalPolygon{T,N}, h :: UInt64) where {N,T <: Integer}
    h = hash(P.rationality, h)
    for i = 1 : N
        h = hash(scaled_vertex(P,i), h)
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


@doc raw"""
    vertices(P :: RationalPolygon)

Return the vertices of a `P`.

"""
vertices(P :: RationalPolygon{T,N}) where {N,T <: Integer} = collect(P)


@doc raw"""
    Base.:(+)(P :: RationalPolygon{T}, Q :: RationalPolygon{T}) where {T <: Integer}

Return the Minkowski sum of two rational polygons sharing the same
rationality.

"""
function Base.:(+)(P :: RationalPolygon{T}, Q :: RationalPolygon{T}) where {T <: Integer}
    rationality(P) == rationality(Q) || error("the rationalities must coincide for the minkowski sum")
    return convex_hull([v + w for v ∈ vertices(P) for w ∈ vertices(Q)], rationality(P))
end


@doc raw"""
    Base.:(+)(P :: RationalPolygon{T}, v :: LatticePoint{T}) where {T <: Integer}

For a `k`-rational polygon `P`, return the translated polygons `P + (v // k)`.

"""
Base.:(+)(P :: RationalPolygon{T}, v :: LatticePoint{T}) where {T <: Integer} =
RationalPolygon(vertex_matrix(P) .+ v, rationality(P))

Base.:(-)(P :: RationalPolygon) =
RationalPolygon(-vertex_matrix(P), rationality(P))

Base.:(-)(P :: RationalPolygon{T}, Q :: RationalPolygon{T}) where {T <: Integer} =
P + (-Q)

Base.:(-)(P :: RationalPolygon{T}, v :: LatticePoint{T}) where {T <: Integer} =
P + (-v)

Base.:(*)(c :: T, P :: RationalPolygon{T}) where {T <: Integer} =
RationalPolygon(c * vertex_matrix(P), rationality(P))

Base.:(*)(c :: Rational{T}, P :: RationalPolygon{T}) where {T <: Integer} =
RationalPolygon(numerator(c) * vertex_matrix(P), denominator(c) * rationality(P))

function Base.:(*)(U :: Matrix2{T}, P :: RationalPolygon{T}) where {T <: Integer}
    det(U) ∈ [1,-1] || error("the given matrix is not unimodular")
    return RationalPolygon(U * vertex_matrix(P), rationality(P))
end

Base.:(//)(P :: RationalPolygon{T}, c :: T) where {T <: Integer} =
RationalPolygon(vertex_matrix(P), c * rationality(P))

Base.:(//)(P :: RationalPolygon{T}, c :: Rational{T}) where {T <: Integer} =
RationalPolygon(denominator(c) * vertex_matrix(P), numerator(c) * rationality(P))

