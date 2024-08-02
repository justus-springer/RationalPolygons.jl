
@doc raw"""
    remove_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

Return a pair `(Q, keeps_genus)`, where `Q` is the convex hull of all
`k`-rational points of `P` except the `i`-th vertex and `keeps_genus :: Bool`
is true if and only if `Q` has the same number of interior lattice points as
`P`. If the argument `primitive = true` is passed, the convex hull is taken of
all _primitive_ `k`-rational points of `P` except the `i`-th vertex.

# Example

Compare Example 10.2.7 of [CLS11](@cite).

```jldoctest
julia> P = convex_hull(LatticePoint{Int}[(0,0),(7,-5),(7,1),(0,1)])
Rational polygon of rationality 1 with 4 vertices.

julia> Q,_ = remove_vertex(P,1)
(Rational polygon of rationality 1 with 4 vertices., false)

julia> vertices(Q)
4-element Vector{SVector{2, Rational{Int64}}}:
 [0, 1]
 [3, -2]
 [7, -5]
 [7, 1]
```

"""
function remove_vertex(P :: RationalPolygon{T,N}, i :: Int; primitive :: Bool = false) where {N,T <: Integer}

    k = rationality(P)
    keeps_genus :: Bool = false

    if primitive
        v0 = vertex(P,i)
        ps = filter(p -> is_primitive(k*p), k_rational_points(P, k))
        filter!(p -> p != v0, ps)
        Q = convex_hull(ps, k)
        keeps_genus = number_of_interior_lattice_points(P) == number_of_interior_lattice_points(Q)
        return (Q, keeps_genus)
    else
        # the shortcut via hilbert bases only works in the non-primitive
        # case at the moment
        u, v, w = scaled_vertex(P, i-1), scaled_vertex(P,i), scaled_vertex(P,i+1)
        p1, p2 = u - v, w - v
        q1, q2 = primitivize(p1), primitivize(p2)
        hb = [p + v for p ∈ hilbert_basis(q1,q2)]
        
        keeps_genus = all(p -> gcd(p) % k != 0, hb[2 : end-1])

        vs = LatticePoint{T}[]
        !is_primitive(p1) && push!(vs, u)
        push!(vs,first(hb))
        for j = 2 : length(hb)-1
            if hb[j]-hb[j-1] != hb[j+1]-hb[j]
                push!(vs, hb[j])
            end
        end
        push!(vs,last(hb))
        !is_primitive(p2) && push!(vs, w)
        for j = i+2 : i+N-2
            push!(vs, scaled_vertex(P,j))
        end

        Q = RationalPolygon(vs, rationality(P))

        return (Q, keeps_genus)

    end

end
