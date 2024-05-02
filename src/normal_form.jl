
@attr function lattice_edge_areas(P :: RationalPolygon)
    n = number_of_vertices(P)
    return [abs(det(P[i+1] - P[i], P[i] - P[i-1])) for i = 1 : n]
end

@attr function special_vertices(P :: RationalPolygon)
    r = number_of_vertices(P)
    ea = lattice_edge_areas(P)
    m = maximum(ea)
    return filter(i -> ea[i] == m, 1 : r)
end

@attr function normal_form(P :: RationalPolygon{T}) where {T <: Integer}

    V = vertex_matrix(P)

    As = Matrix{T}[]
    for i ∈ special_vertices(P)
        push!(As, [V[:,i:end] V[:,begin:i-1]])
        push!(As, [V[:,i:-1:begin] V[:,end:-1:i+1]])
    end
    As = map(hnf, As)

    # take the lexicographical minumum of all the hermite normal forms
    _lt(A, B) = vec(A) < vec(B)
    A = sort(As; lt = _lt)[1]

    Q = RationalPolygon(A, rationality(P)) 

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
