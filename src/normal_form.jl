
@attr function lattice_edges(P :: RationalPolygon)
    vs, r = lattice_vertices(P), number_of_vertices(P)
    return [(vs[i], vs[mod(i+1,1:r)]) for i = 1 : r]
end

@attr lattice_edge_vectors(P :: RationalPolygon) =
[e[2] - e[1] for e ∈ lattice_edges(P)]

@attr function lattice_edge_areas(P :: RationalPolygon)
    ev = lattice_edge_vectors(P)
    r = number_of_vertices(P)
    return [abs(det(ev[mod(i-1,1:r)], ev[i])) for i = 1 : r]
end

@attr function special_vertices(P :: RationalPolygon)
    r = number_of_vertices(P)
    ea = lattice_edge_areas(P)
    m = maximum(ea)
    return filter(i -> ea[i] == m, 1 : r)
end

@attr function normal_form(P :: RationalPolygon{T}) where {T <: Integer}
    n = length(vertices(P))
    vs = lattice_vertices(P)
    k = rationality(P)

    cycls = [map(x -> mod(i+x, 1:n), 0:n-1) for i ∈ special_vertices(P)]
    append!(cycls, [map(x -> mod(i-x, 1:n), 0:n-1) for i ∈ special_vertices(P)])

    # take the hermite normal form of all possible permutations of the
    # rays keeping the counterclockwise ordering
    As = map(is -> hcat(map(v -> [v[1],v[2]], vs[is])...), cycls)
    As = map(hnf, As)
    # take the lexicographical minumum of all the hermite normal forms
    _lt(A, B) = vcat(A...) < vcat(B...)
    A = sort(As; lt = _lt)[1]

    Q = ConvexHull([(T(A[1,i])//k,T(A[2,i])//k) for i = 1 : n]; rationality = rationality(P))

    ### attribute carrying ###
    set_attribute!(Q, :is_normal_form, true)
    set_attribute!(Q, :normal_form, Q)
    if has_attribute(P, :number_of_interior_lattice_points)
        set_attribute!(Q, :number_of_interior_lattice_points, number_of_interior_lattice_points(P))
    end

    return Q

end

@attr is_normal_form(P :: RationalPolygon) =
P == normal_form(P)

function are_equivalent(P :: RationalPolygon, Q :: RationalPolygon)
    return normal_form(P) == normal_form(Q)
end
