
@doc raw"""
    classify_maximal_polygons_genus_one_m1p1!(k :: T; Ps :: Vector{RationalPolygon{T}} = RationalPolygon{T}[]) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point that can be realized in $\mathbb{Q} \times [-1,1]$.
Optionally, a set `Ps` of starting polygons can be given, in which case
the resulting polygons will be added to `Ps`, if they are not already
present.

"""
function classify_maximal_polygons_genus_one_m1p1!(k :: T; Ps :: Vector{RationalPolygon{T}} = RationalPolygon{T}[], primitive :: Bool = false) where {T <: Integer}
    A = convex_hull([(-k,zero(T)),(zero(T),k),(-k,k)], k)
    B = convex_hull([(-k,zero(T)),(k,zero(T)),(2k,k),(-k,k)], k)

    vs = [v for v ∈ k_rational_points(k,A; primitive) if v[2] > 0]
    ws = [w for w ∈ k_rational_points(k,B; primitive) if w[2] > 0]

    for v ∈ vs, w ∈ ws
        H1 = affine_halfplane(v,(-T(1),T(0)))
        H2 = affine_halfplane((T(1),T(0)),w)
        w ∈ H1 && v ∈ H2 || continue

        H_upper = affine_halfplane((T(0),-T(1)),-T(1))
        H_lower = affine_halfplane((T(0),T(1)),-T(1))

        P = k_rational_hull(k, intersect_halfplanes([H1,H2,H_upper,H_lower]); primitive)
        interior_lattice_points(P) == [(0,0)] || continue
        all(Q -> !are_equivalent(P,Q), Ps) || continue
        is_maximal(P; primitive) || continue

        push!(Ps, normal_form(P))

    end

    return Ps
end

classify_maximal_polygons_genus_one_m1p1(k :: T; primitive :: Bool = false) where {T <: Integer} =
classify_maximal_polygons_genus_one_m1p1!(k; Ps = RationalPolygon{T}[], primitive)


@doc raw"""
    classify_maximal_polygons_genus_one_m1p2!(k :: T; Ps :: Vector{RationalPolygon{T}} = RationalPolygon{T}[]) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point that can be realized in $\mathbb{Q} \times [-1,2]$.
Optionally, a set `Ps` of starting polygons can be given, in which case
the resulting polygons will be added to `Ps`, if they are not already
present.

"""
function classify_maximal_polygons_genus_one_m1p2!(k :: T; Ps :: Vector{RationalPolygon{T}} = RationalPolygon{T}[], primitive :: Bool = false) where {T <: Integer}
    A = convex_hull([(-k,k), (zero(T),k), (zero(T),2k), (-2k, 2k)], k)
    B = convex_hull([(-k,zero(T)), (k,zero(T)), (3k,-k), (-2k,-k)], k)

    a1,a2 = (-T(1),T(1)), (T(0),T(1))
    b1,b2 = (-T(1),T(0)), (T(1),T(0))
    
    vs = filter(v -> v[2] > 1, k_rational_points(k, A; primitive))
    ws = filter(w -> w[2] < 0, k_rational_points(k, B; primitive))

    for v1 ∈ vs, v2 ∈ vs
        Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
        v1 ∈ Ha2 && v2 ∈ Ha1 || continue

        for w1 ∈ ws, w2 ∈ ws
            w1 ∈ Ha1 && w1 ∈ Ha2 && w2 ∈ Ha1 && w2 ∈ Ha2 || continue

            Hb1, Hb2 = affine_halfplane(b1,w1), affine_halfplane(w2,b2)
            w1 ∈ Hb2 && w2 ∈ Hb1 || continue
            v1 ∈ Hb1 && v1 ∈ Hb2 && v2 ∈ Hb1 && v2 ∈ Hb2 || continue

            H_upper = affine_halfplane((T(1),T(2)),(T(0),T(2)))
            H_lower = affine_halfplane((T(0),-T(1)),(T(1),-T(1)))

            P = k_rational_hull(k, intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]); primitive)
            interior_lattice_points(P) == [(0,0)] || continue
            all(Q -> !are_equivalent(P,Q), Ps) || continue
            is_maximal(P; primitive) || continue

            push!(Ps, normal_form(P))
        end
    end

    return Ps

end

classify_maximal_polygons_genus_one_m1p2(k :: T; primitive :: Bool = false) where {T <: Integer} =
classify_maximal_polygons_genus_one_m1p2!(k; Ps = RationalPolygon{T}[], primitive)

@doc raw"""
    classify_maximal_polygons_genus_one(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point.

"""
function classify_maximal_polygons_genus_one(k :: T ; primitive :: Bool = false, logging :: Bool = false) where {T <: Integer}

    primstring = primitive ? "primitive " : ""

    Ps = RationalPolygon{T}[]
    
    Ps = classify_maximal_polygons_genus_one_m1p1!(k; Ps, primitive)
    total_count = length(Ps)

    logging && @info "Found $(length(Ps)) $(primstring)maximal polygons in QQ x [-1,1]."
    total_count = length(Ps)

    Ps = classify_maximal_polygons_genus_one_m1p2!(k; Ps, primitive)

    logging && @info "Found $(length(Ps) - total_count) new $(primstring)maximal polygons in QQ x [-1,2]. Total : $(length(Ps))"
    total_count = length(Ps)

    logging && @info "(skipping [-2,2])"
    logging && @info "Found $(length(Ps) - total_count) new $(primstring)maximal polygons in QQ x [-2,2]. Total : $(length(Ps))"

    return Ps
end

function classify_polygons_genus_one(k :: T; primitive :: Bool = false, out_path :: Union{Missing,String} = missing, logging = false) where {T <: Integer}

    primstring = primitive ? "primitive " : ""
    logging && @info "Beginning classification of all $(primstring)$k-rational polygons with one interior lattice point."

    Ps = classify_maximal_polygons_genus_one(k; primitive, logging)
    return subpolygons(Ps; primitive, out_path, logging)
end

