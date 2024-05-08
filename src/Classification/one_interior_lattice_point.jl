
@doc raw"""
    classify_maximal_polygons_genus_one_m1p1(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point that can be realized in $\mathbb{Q} \times [-1,1]$.

"""
function classify_maximal_polygons_genus_one_m1p1(k :: T) where {T <: Integer}
    Ps = RationalPolygon{T}[]

    A = convex_hull([(-k,zero(T)),(zero(T),k),(-k,k)], k)
    B = convex_hull([(-k,zero(T)),(k,zero(T)),(2k,k),(-k,k)], k)

    vs = [v for v ∈ k_rational_points(k,A) if v[2] > 0]
    ws = [w for w ∈ k_rational_points(k,B) if w[2] > 0]

    for v ∈ vs, w ∈ ws
        H1 = affine_halfplane(v,(-T(1),T(0)))
        H2 = affine_halfplane((T(1),T(0)),w)
        w ∈ H1 && v ∈ H2 || continue

        H_upper = affine_halfplane((T(0),-T(1)),-T(1))
        H_lower = affine_halfplane((T(0),T(1)),-T(1))

        P = k_rational_hull(k, intersect_halfplanes([H1,H2,H_upper,H_lower]))
        interior_lattice_points(P) == [(0,0)] || continue
        all(Q -> !are_equivalent(P,Q), Ps) || continue
        is_maximal(P) || continue

        push!(Ps, normal_form(P))

    end

    return Ps
end


@doc raw"""
    classify_maximal_polygons_genus_one_m1p2(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point that can be realized in $\mathbb{Q} \times [-1,2]$.

"""
function classify_maximal_polygons_genus_one_m1p2(k :: T) where {T <: Integer}
    Ps = RationalPolygon{T}[]

    A = convex_hull([(-k,k), (zero(T),k), (zero(T),2k), (-2k, 2k)], k)
    B = convex_hull([(-k,zero(T)), (k,zero(T)), (3k,-k), (-2k,-k)], k)

    a1,a2 = (-T(1),T(1)), (T(0),T(1))
    b1,b2 = (-T(1),T(0)), (T(1),T(0))
    
    vs = filter(v -> v[2] > 1, k_rational_points(k, A))
    ws = filter(w -> w[2] < 0, k_rational_points(k, B))

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

            P = k_rational_hull(k, intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]))
            interior_lattice_points(P) == [(0,0)] || continue
            all(Q -> !are_equivalent(P,Q), Ps) || continue
            is_maximal(P) || continue

            push!(Ps, normal_form(P))
        end
    end

    return Ps

end

@doc raw"""
    classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{Q} \times [-2,2]$ that have non-empty
intersection with the `q`-th classification box, where `1 ≤ q ≤ 5`.

"""
function classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int) where {T <: Integer}
    
    Ps = RationalPolygon{T}[]

    A = convex_hull([(-k,k), (zero(T),k), (zero(T),2k), (-2k, 2k)], k)
    B = convex_hull([(-k,zero(T)), (k,zero(T)), (3k,-k), (-2k,-k)], k)

    a1,a2 = (-one(T),one(T)), (zero(T),one(T))
    b1,b2 = (-one(T),zero(T)), (one(T),zero(T))

    if q == 1
        C = convex_hull([(-2k,-k),(-3k,-2k),(-2k,-2k),(-k,-k)], k)
        c1,c2 = (-2*one(T), -one(T)), (-one(T), -one(T))
    elseif q == 2
        C = convex_hull([(-k,-k),(-2k,-2k),(T(0),-2k),(T(0),-k)], k)
        c1,c2 = (-one(T), -one(T)), (zero(T), -one(T))
    elseif q == 3
        C = convex_hull([(T(0),-k),(T(0),-2k),(2k,-2k),(k,-k)], k)
        c1,c2 = (zero(T), -one(T)), (one(T), -one(T))
    elseif q == 4
        C = convex_hull([(k,-k),(2k,-2k),(4k,-2k),(2k,-k)], k)
        c1,c2 = (one(T), -one(T)), (2*one(T), -one(T))
    elseif q == 5
        C = convex_hull([(2k,-k),(4k,-2k),(5k,-2k),(3k,-k)], k)
        c1,c2 = (2*one(T), -one(T)), (3*one(T), -one(T))
    end
    
    vs = filter(v -> v[2] > 1, k_rational_points(k, A))
    ws = unique(filter(w -> w[2] < 0, [k_rational_points(k, B) ; k_rational_points(k,C)]))
    us = filter(u -> u[2] < -1, k_rational_points(k,C))

    for v1 ∈ vs, v2 ∈ vs
        Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
        v1 ∈ Ha2 && v2 ∈ Ha1 || continue
        
        for w1 ∈ ws, w2 ∈ ws
            w1 ∈ Ha1 && w1 ∈ Ha2 && w2 ∈ Ha1 && w2 ∈ Ha2 || continue

            Hb1, Hb2 = affine_halfplane(b1,w1), affine_halfplane(w2,b2)
            w1 ∈ Hb2 && w2 ∈ Hb1 || continue
            v1 ∈ Hb1 && v1 ∈ Hb2 && v2 ∈ Hb1 && v2 ∈ Hb2 || continue

            for u1 ∈ us, u2 ∈ us

                u1 ∈ Ha1 && u1 ∈ Ha2 && u2 ∈ Ha1 && u2 ∈ Ha2 || continue
                u1 ∈ Hb1 && u1 ∈ Hb2 && u2 ∈ Hb1 && u2 ∈ Hb2 || continue

                Hc1, Hc2 = affine_halfplane(c1,u1), affine_halfplane(u2,c2)
                u1 ∈ Hc2 && u2 ∈ Hc1 || continue

                v1 ∈ Hc1 && v1 ∈ Hc2 && v2 ∈ Hc1 && v2 ∈ Hc2 || continue
                w1 ∈ Hc1 && w1 ∈ Hc2 && w2 ∈ Hc1 && w2 ∈ Hc2 || continue

                H_upper = affine_halfplane((T(1),T(2)),(T(0),T(2)))
                H_lower = affine_halfplane((T(0),-T(2)),(T(1),-T(2)))

                P = k_rational_hull(k, intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,Hc1,Hc2,H_upper,H_lower]))
                interior_lattice_points(P) == [(0,0)] || continue
                all(Q -> !are_equivalent(P,Q), Ps) || continue
                is_maximal(P) || continue

                push!(Ps, normal_form(P))
            end
        end
    end

    return Ps

end


@doc raw"""
    classify_maximal_polygons_genus_one(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point.

"""
function classify_maximal_polygons_genus_one(k :: T ; logging = false) where {T <: Integer}
    Ps = RationalPolygon{T}[]
    count, total_count = 0, 0
    
    new_Ps = classify_maximal_polygons_genus_one_m1p1(k)
    count = length(new_Ps)
    unique!(append!(Ps, new_Ps))
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count maximal polygons in QQ x [-1,1]. New: $new_count, total: $total_count"

    new_Ps = classify_maximal_polygons_genus_one_m1p2(k)
    count = length(new_Ps)
    unique!(append!(Ps, new_Ps))
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count maximal polygons in QQ x [-1,2]. New: $new_count, total: $total_count"
    
    for q = 1 : 5

        new_Ps = classify_maximal_polygons_genus_one_m2p2(k, q)
        count = length(new_Ps)
        unique!(append!(Ps, new_Ps))
        new_count = length(Ps) - total_count
        total_count = length(Ps)

        logging && @info "Found $count maximal polygons in QQ x [-2,2], box $q. New: $new_count, total: $total_count"

    end

    return Ps
end

function classify_polygons_genus_one(k :: T; out_path :: Union{Missing,String} = missing, logging = false) where {T <: Integer}

    logging && @info "Beginning classification of all $k-rational polygons with one interior lattice point."

    Ps = classify_maximal_polygons_genus_one(k; logging)
    return subpolygons(Ps; out_path, logging)
end

