
@doc raw"""
    classify_maximal_polygons_genus_one_m1p1(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point that can be realized in $\mathbb{Q} \times [-1,1]$.

"""
function classify_maximal_polygons_genus_one_m1p1(k :: T; primitive :: Bool = false) where {T <: Integer}
    Ps = RationalPolygon{T}[]

    A = convex_hull(RationalPoint{T}.([(-1,0),(0,1),(-1,1)]), k)
    B = convex_hull(RationalPoint{T}.([(-1,0),(1,0),(2,1),(-1,1)]), k)

    vs = [v for v ∈ k_rational_points(k,A; primitive) if v[2] > 0]
    ws = [w for w ∈ k_rational_points(k,B; primitive) if w[2] > 0]

    for v ∈ vs, w ∈ ws
        H1 = affine_halfplane(v,RationalPoint{T}(-1,0))
        H2 = affine_halfplane(RationalPoint{T}(1,0),w)
        w ∈ H1 && v ∈ H2 || continue

        H_upper = affine_halfplane(RationalPoint{T}(0,-1),-T(1))
        H_lower = affine_halfplane(RationalPoint{T}(0,1),-T(1))

        P = k_rational_hull(k, intersect_halfplanes([H1,H2,H_upper,H_lower]); primitive)
        interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
        all(Q -> !are_equivalent(P,Q), Ps) || continue

        !primitive && (is_maximal(P) || continue)

        push!(Ps, normal_form(P))

    end

    return Ps
end


@doc raw"""
    classify_maximal_polygons_genus_one_m1p2(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior
lattice point that can be realized in $\mathbb{Q} \times [-1,2]$.

"""
function classify_maximal_polygons_genus_one_m1p2(k :: T; primitive :: Bool = false) where {T <: Integer}
    Ps = RationalPolygon{T}[]

    A = convex_hull(RationalPoint{T}.([(-1,1),(0,1),(0,2),(-2,2)]), k)
    B = convex_hull(RationalPoint{T}.([(-1,0),(1,0),(3,-1),(-2,-1)]), k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(1,0)
    
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

            H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
            H_lower = affine_halfplane(RationalPoint{T}(0,-1),RationalPoint{T}(1,-1))

            P = k_rational_hull(k, intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]); primitive)
            interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
            all(Q -> !are_equivalent(P,Q), Ps) || continue

            !primitive && (is_maximal(P) || continue)

            push!(Ps, normal_form(P))
        end
    end

    return Ps

end

@doc raw"""
    classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{Q} \times [-2,2]$ that have non-empty
intersection with the `q`-th classification box, where `1 ≤ q ≤ 3`.

"""
function classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int; primitive :: Bool = false) where {T <: Integer}
    
    Ps = RationalPolygon{T}[]

    A = convex_hull(RationalPoint{T}[(-1,1),(0,1),(0,2),(-2,2)], k)
    B = convex_hull(RationalPoint{T}[(-1,0),(1,0),(3,-1),(-2,-1)], k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(1,0)

    if q == 1
        B1 = convex_hull(RationalPoint{T}[(-2,-1),(-1,-1),(-1,0)], k)
        B2 = convex_hull(RationalPoint{T}[(-2,-1),(3,-1),(1,0)], k)
        C = convex_hull(RationalPoint{T}[(-2,-1),(-3,-2),(-2,-2),(-1,-1)], k)
        c1,c2 = RationalPoint{T}(-2,-1), RationalPoint{T}(-1,-1)
    elseif q == 2
        B1 = convex_hull(RationalPoint{T}[(-2,-1),(0,-1),(-1,0)], k)
        B2 = convex_hull(RationalPoint{T}[(-1,-1),(3,-1),(1,0)], k)
        C = convex_hull(RationalPoint{T}[(-1,-1),(-2,-2),(0,-2),(0,-1)], k)
        c1,c2 = (-one(T), -one(T)), (zero(T), -one(T))
        c1,c2 = RationalPoint{T}(-1,-1), RationalPoint{T}(0,-1)
    elseif q == 3
        B1 = convex_hull(RationalPoint{T}[(-2,-1),(1,-1),(-1,0)], k)
        B2 = convex_hull(RationalPoint{T}[(0,-1),(3,-1),(1,0)], k)
        C = convex_hull(RationalPoint{T}[(0,-1),(0,-2),(2,-2),(1,-1)], k)
        c1,c2 = RationalPoint{T}(0,-1), RationalPoint{T}(1,-1)
    end
    
    vs = k_rational_points(k,A; primitive)
    filter!(v -> v[2] > 1, vs)

    us = k_rational_points(k,C; primitive)
    filter!(u -> u[2] < -1, us)

    ws1 = append!(k_rational_points(k,B1; primitive), us)
    filter!(w -> w[2] < 0, ws1)

    ws2 = append!(k_rational_points(k,B2; primitive), us)
    filter!(w -> w[2] < 0, ws2)

    for v1 ∈ vs, v2 ∈ vs
        Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
        v1 ∈ Ha2 && v2 ∈ Ha1 || continue
        
        for w1 ∈ ws1, w2 ∈ ws2
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

                H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
                H_lower = affine_halfplane(RationalPoint{T}(0,-2),RationalPoint{T}(1,-2))

                P = k_rational_hull(k, intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,Hc1,Hc2,H_upper,H_lower]); primitive)
                interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
                all(Q -> !are_equivalent(P,Q), Ps) || continue

                !primitive && (is_maximal(P) || continue)

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
function classify_maximal_polygons_genus_one(k :: T ; primitive :: Bool = false, logging :: Bool = false) where {T <: Integer}

    primstring = primitive ? "primitive " : ""

    Ps = RationalPolygon{T}[]
    count, total_count = 0, 0
    
    new_Ps = classify_maximal_polygons_genus_one_m1p1(k; primitive)
    count = length(new_Ps)
    unique!(append!(Ps, new_Ps))
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count $(primstring)maximal polygons in QQ x [-1,1]. New: $new_count, total: $total_count"

    new_Ps = classify_maximal_polygons_genus_one_m1p2(k; primitive)
    count = length(new_Ps)
    unique!(append!(Ps, new_Ps))
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count $(primstring)maximal polygons in QQ x [-1,2]. New: $new_count, total: $total_count"
    
    for q = 1 : 3

        new_Ps = classify_maximal_polygons_genus_one_m2p2(k, q; primitive)
        count = length(new_Ps)
        unique!(append!(Ps, new_Ps))
        new_count = length(Ps) - total_count
        total_count = length(Ps)

        logging && @info "Found $count $(primstring)maximal polygons in QQ x [-2,2], box $q. New: $new_count, total: $total_count"

    end

    return Ps
end

function classify_polygons_genus_one(k :: T; primitive :: Bool = false, out_path :: Union{Missing,String} = missing, logging = false) where {T <: Integer}

    primstring = primitive ? "primitive " : ""
    logging && @info "Beginning classification of all $(primstring)$k-rational polygons with one interior lattice point."

    Ps = classify_maximal_polygons_genus_one(k; primitive, logging)
    return subpolygons(Ps; primitive, out_path, logging)
end

@doc raw"""
    classify_maximal_lattice_free_polygons_m1p1(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with no interior lattice points that
can be realized in $\mathbb{Q} \times [-1,1]$.

"""
function classify_maximal_lattice_free_polygons_m1p1(k :: T) where {T <: Integer}
    Ps = RationalPolygon{T}[]

    A = convex_hull(RationalPoint{T}.([(-1,0),(0,1),(-1,1)]), k)
    B = convex_hull(RationalPoint{T}.([(-1,0),(0,0),(1,1),(-1,1)]), k)

    vs = [v for v ∈ k_rational_points(k,A) if v[2] > 0]
    ws = [w for w ∈ k_rational_points(k,B) if w[2] > 0]

    for v ∈ vs, w ∈ ws
        H1 = affine_halfplane(v,RationalPoint{T}(-1,0))
        H2 = affine_halfplane(RationalPoint{T}(0,0),w)
        w ∈ H1 && v ∈ H2 || continue

        H_upper = affine_halfplane(RationalPoint{T}(0,-1),-T(1))
        H_lower = affine_halfplane(RationalPoint{T}(0,1),-T(1))

        P = k_rational_hull(k, intersect_halfplanes([H1,H2,H_upper,H_lower]))
        number_of_interior_lattice_points(P) == 0 || continue
        all(Q -> !are_equivalent(P,Q), Ps) || continue
        is_maximal(P) || continue

        push!(Ps, normal_form(P))

    end

    return Ps
    
end


@doc raw"""
    classify_maximal_lattice_free_polygons_m1p2(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with no interior lattice points that
can be realized in $\mathbb{Q} \times [-1,2]$.

"""
function classify_maximal_lattice_free_polygons_m1p2(k :: T) where {T <: Integer}
    Ps = RationalPolygon{T}[]

    A = convex_hull(RationalPoint{T}.([(0,1),(1,2),(-2,2),(-1,1)]), k)
    B = convex_hull(RationalPoint{T}.([(0,0),(1,-1),(-2,-1),(-1,0)]), k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(0,0)
    
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

            H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
            H_lower = affine_halfplane(RationalPoint{T}(0,-1),RationalPoint{T}(1,-1))

            P = k_rational_hull(k, intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]))
            number_of_vertices(P) > 2 || continue
            number_of_interior_lattice_points(P) == 0 || continue
            all(Q -> !are_equivalent(P,Q), Ps) || continue
            is_maximal(P) || continue

            push!(Ps, normal_form(P))
        end
    end

    return Ps
end


@doc raw"""
    classify_maximal_lattice_free_polygons(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with no interior lattice points.

"""
function classify_maximal_lattice_free_polygons(k :: T ; logging = false) where {T <: Integer}
    Ps = RationalPolygon{T}[]
    count, total_count = 0, 0
    
    new_Ps = classify_maximal_lattice_free_polygons_m1p1(k)
    count = length(new_Ps)
    unique!(append!(Ps, new_Ps))
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count lattice-free maximal polygons in QQ x [-1,1]. New: $new_count, total: $total_count"

    new_Ps = classify_maximal_lattice_free_polygons_m1p2(k)
    count = length(new_Ps)
    unique!(append!(Ps, new_Ps))
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count lattice-free maximal polygons in QQ x [-1,2]. New: $new_count, total: $total_count"

    return Ps
end