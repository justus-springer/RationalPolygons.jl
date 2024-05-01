@doc raw"""
    remove_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

Remove the `i`-th vertex from `P`, i.e. return the convex hull of all
`k`-rational points of `P` except the `i`-th vertex.

"""
function remove_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}
    k = rationality(P)
    vs, r = lattice_vertices(P), number_of_vertices(P)
    p1 = vs[mod(i-1,1:r)] - vs[i]
    q1 = primitivize(p1)
    p2 = vs[mod(i+1,1:r)] - vs[i]
    q2 = primitivize(p2)
    hb = [v + vs[i] for v ∈ hilbert_basis(q1,q2)]


    new_vs = empty(vs)
    !is_primitive(p1) && push!(new_vs, vs[mod(i-1,1:r)])
    push!(new_vs,first(hb))
    for j = 2 : length(hb)-1
        if hb[j]-hb[j-1] != hb[j+1]-hb[j]
            push!(new_vs, hb[j])
        end
    end
    push!(new_vs,last(hb))
    !is_primitive(p2) && push!(new_vs, vs[mod(i+1,1:r)])
    for j = i+2 : i+r-2
        push!(new_vs, vs[mod(j,1:r)])
    end

    Q = ConvexHull(new_vs .// k; rationality = k)

    ### attribute carrying ###
    
    if has_attribute(P, :interior_lattice_points)
        new_int_points = filter(p -> p ∉ hb, interior_lattice_points(P))
        set_attribute!(Q, :interior_lattice_points, new_int_points)
        set_attribute!(Q, :number_of_interior_lattice_points, length(new_int_points))
    end

    return Q

end

@doc raw"""
    subpolygons!(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

Given some rational polygons `Ps` with shared rationality `k` and
number of interior lattice points `n`, add all subpolygons to `Ps`
with the same rationality and number of interior lattice points.

"""
function subpolygons!(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    k = rationality(first(Ps))
    n = number_of_interior_lattice_points(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("all polygons must have the same rationality")
    all(P -> number_of_interior_lattice_points(P) == n, Ps) || error("all polygons must have the same number of interior lattice points")

    vertices_dict = Dict{Int,Set{RationalPolygon{T}}}()
    i = 1
    while i ≤ length(Ps)
        P = Ps[i]
        for j = 1 : number_of_vertices(P)
            P2 = normal_form(remove_vertex(P, j))
            number_of_vertices(P2) > 2 || continue
            number_of_interior_lattice_points(P2) == n || continue

            nv = number_of_vertices(P2)
            if !haskey(vertices_dict, nv)
                vertices_dict[nv] = Set{RationalPolygon{T}}()
            end
            P2 ∉ vertices_dict[nv] || continue
            push!(vertices_dict[nv], P2)

            push!(Ps, P2)

        end

        i += 1

    end

    return Ps
                    
end

subpolygons(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer} =
subpolygons!(deepcopy(Ps))

subpolygons(P :: RationalPolygon{T}) where {T <: Integer} =
subpolygons!([P])
