
_cycles(n :: Int) = [map(x -> mod(x+i, 1:n), 1:n) for i = 1:n]

@attr function normal_form(P :: RationalPolygon{T}) where {T <: Integer}
    n = length(vertices(P))
    vs = lattice_vertices(P)
    k = rationality(P)

    # take the hermite normal form of all possible permutations of the
    # rays keeping the counterclockwise ordering
    cycls = [_cycles(n) ; map(reverse, _cycles(n))]
    As = map(is -> hnf(matrix(ZZ, hcat(map(v -> [v[1],v[2]], vs[is])...))), cycls)
    # take the lexicographical minumum of all the hermite normal forms
    _lt(A, B) = vcat(A...) < vcat(B...)
    A = sort(As; lt = _lt)[1]

    Q = ConvexHull([(T(A[1,i])//k,T(A[2,i])//k) for i = 1 : n]; rationality = rationality(P))

    set_attribute!(Q, :is_normal_form, true)
    set_attribute!(Q, :normal_form, Q)

    return Q
end

@attr is_normal_form(P :: RationalPolygon) =
P == normal_form(P)

are_equivalent(P :: RationalPolygon, Q :: RationalPolygon) =
normal_form(P) == normal_form(Q)
