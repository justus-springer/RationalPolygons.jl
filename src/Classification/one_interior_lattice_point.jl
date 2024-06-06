
@doc raw"""
    classify_maximal_polygons_genus_one_m1p1(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{Q} \times [-1,1]$. If `primtive = true`
is passed, then only primitive polygons (i.e. fano polygons) are returned.

"""
function classify_maximal_polygons_genus_one_m1p1(k :: T; primitive :: Bool = false) where {T <: Integer}

    A = RationalPolygon(SA[-k 0 -k ; 0 k k], k)
    B = RationalPolygon(SA[-k k 2k -k ; 0 0 k k], k)

    vs = filter(v -> v[2] > 0, k_rational_points(A,k))
    primitive && filter!(v -> is_primitive(k*v), vs)
    ws = filter(w -> w[2] > 0, k_rational_points(B,k))
    primitive && filter!(w -> is_primitive(k*w), ws)

    Pss = Vector{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, RationalPolygon{T}[])
    end

    Threads.@threads for v ∈ vs
        for w ∈ ws
            H1 = affine_halfplane(v,RationalPoint{T}(-1,0))
            H2 = affine_halfplane(RationalPoint{T}(1,0),w)
            w ∈ H1 && v ∈ H2 || continue

            H_upper = affine_halfplane(RationalPoint{T}(0,-1),-T(1))
            H_lower = affine_halfplane(RationalPoint{T}(0,1),-T(1))

            P = k_rational_hull(intersect_halfplanes([H1,H2,H_upper,H_lower]), k; primitive)
            interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
            all(Q -> !are_unimodular_equivalent(P,Q), Pss[Threads.threadid()]) || continue

            !primitive && (is_maximal(P) || continue)

            push!(Pss[Threads.threadid()], unimodular_normal_form(P))
        end
    end

    Ps = RationalPolygon{T}[]
    for k = 1 : Threads.nthreads()
        for P ∈ Pss[k]
            all(Q -> !are_unimodular_equivalent(P,Q), Ps) || continue
            push!(Ps, P)
        end
    end

    return Ps
end


@doc raw"""
    classify_maximal_polygons_genus_one_m1p2(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{Q} \times [-1,2]$. If `primtive = true`
is passed, then only primitive polygons (i.e. fano polygons) are returned.

"""
function classify_maximal_polygons_genus_one_m1p2(k :: T; primitive :: Bool = false) where {T <: Integer}

    A = RationalPolygon(SA[-2k -k 0 0 ; 2k k k 2k], k)
    B = RationalPolygon(SA[-2k 3k k -k ; -k -k 0 0], k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(1,0)
    
    vs = filter(v -> v[2] > 1, k_rational_points(A, k))
    primitive && filter!(v -> is_primitive(k*v), vs)
    ws = filter(w -> w[2] < 0, k_rational_points(B, k))
    primitive && filter!(w -> is_primitive(k*w), ws)

    Pss = Vector{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, RationalPolygon{T}[])
    end

    Threads.@threads for v1 ∈ vs
        for v2 ∈ vs
            Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
            v1 ∈ Ha2 && v2 ∈ Ha1 || continue

            for w1 ∈ ws, w2 ∈ ws
                w1 ∈ Ha1 && w1 ∈ Ha2 && w2 ∈ Ha1 && w2 ∈ Ha2 || continue

                Hb1, Hb2 = affine_halfplane(b1,w1), affine_halfplane(w2,b2)
                w1 ∈ Hb2 && w2 ∈ Hb1 || continue
                v1 ∈ Hb1 && v1 ∈ Hb2 && v2 ∈ Hb1 && v2 ∈ Hb2 || continue

                H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
                H_lower = affine_halfplane(RationalPoint{T}(0,-1),RationalPoint{T}(1,-1))

                P = k_rational_hull(intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]), k; primitive)
                interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
                all(Q -> !are_unimodular_equivalent(P,Q), Pss[Threads.threadid()]) || continue

                !primitive && (is_maximal(P) || continue)

                push!(Pss[Threads.threadid()], unimodular_normal_form(P))
            end
        end
    end

    Ps = RationalPolygon{T}[]
    for k = 1 : Threads.nthreads()
        for P ∈ Pss[k]
            all(Q -> !are_unimodular_equivalent(P,Q), Ps) || continue
            push!(Ps, P)
        end
    end

    return Ps

end

@doc raw"""
    classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{Q} \times [-2,2]$ that have non-empty
intersection with the `q`-th classification box, where `1 ≤ q ≤ 3`. If
`primtive = true` is passed, then only primitive polygons (i.e. fano polygons)
are returned.

"""
function classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int; primitive :: Bool = false) where {T <: Integer}
    

    A = RationalPolygon(SA[-2k -k 0 0 ; 2k k k 2k], k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(1,0)

    if q == 1
        B1 = RationalPolygon(SA[-2k -k -k ; -k -k 0], k)
        B2 = RationalPolygon(SA[-2k 3k k ; -k -k 0], k)
        C = RationalPolygon(SA[-3k -2k -k -2k ; -2k -2k -k -k], k)
        c1,c2 = RationalPoint{T}(-2,-1), RationalPoint{T}(-1,-1)
    elseif q == 2
        B1 = RationalPolygon(SA[-2k 0 -k ; -k -k 0], k)
        B2 = RationalPolygon(SA[-k 3k k ; -k -k 0], k)
        C = RationalPolygon(SA[-2k 0 0 -k ; -2k -2k -k -k], k)
        c1,c2 = (-one(T), -one(T)), (zero(T), -one(T))
        c1,c2 = RationalPoint{T}(-1,-1), RationalPoint{T}(0,-1)
    elseif q == 3
        B1 = RationalPolygon(SA[-2k k -k ; -k -k 0], k)
        B2 = RationalPolygon(SA[0 3k k ; -k -k 0], k)
        C = RationalPolygon(SA[0 2k k 0 ; -2k -2k -k -k], k)
        c1,c2 = RationalPoint{T}(0,-1), RationalPoint{T}(1,-1)
    end
    
    vs = filter(v -> v[2] > 1, k_rational_points(A,k))
    primitive && filter!(v -> is_primitive(k*v), vs)

    us = filter(u -> u[2] < -1, k_rational_points(C,k))
    primitive && filter!(u -> is_primitive(k*u), us)

    ws1 = append!(filter(w -> w[2] < 0, k_rational_points(B1,k)), us)
    primitive && filter!(w -> is_primitive(k*w), ws1)

    ws2 = append!(filter(w -> w[2] < 0, k_rational_points(B2,k)), us)
    primitive && filter!(w -> is_primitive(k*w), ws2)

    Pss = Vector{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, RationalPolygon{T}[])
    end

    Threads.@threads for v1 ∈ vs
        for v2 ∈ vs
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

                    P = k_rational_hull(intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,Hc1,Hc2,H_upper,H_lower]), k; primitive)
                    interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
                    all(Q -> !are_unimodular_equivalent(P,Q), Pss[Threads.threadid()]) || continue

                    !primitive && (is_maximal(P) || continue)

                    push!(Pss[Threads.threadid()], unimodular_normal_form(P))
                end
            end
        end
    end

    Ps = RationalPolygon{T}[]
    for k = 1 : Threads.nthreads()
        for P ∈ Pss[k]
            all(Q -> !are_unimodular_equivalent(P,Q), Ps) || continue
            push!(Ps, P)
        end
    end

    return Ps

end


@doc raw"""
    classify_maximal_polygons_genus_one(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point. If `primtive = true` is passed, then only primitive polygons (i.e. fano
polygons) are returned.

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


@doc raw"""
    classify_polygons_genus_one(k :: T) where {T <: Integer}

Compute all `k`-rational polygons with exactly one interior lattice point. The following keyword arguments are supported:

- `primitive :: Bool`. If set to true, only primitive polygons (i.e. fano polygons) are returned.

- `out_path :: Union{Missing,String}`. If set to `missing`, all polygons will
be kept in memory and returned. If set to a path to an empty directory, the
storage of the polygons will be delegated to the disk, where the polygons will
be saved into text files according to their normalized area. Since the number
of polygons grows rapidly in `k`, this is the recommended method for `k ≥ 4`.

- `logging :: Bool`. Controls whether to display logging messages showing the
current progress.

"""
function classify_polygons_genus_one(k :: T; primitive :: Bool = false, out_path :: Union{Missing,String} = missing, logging = false) where {T <: Integer}

    primstring = primitive ? "primitive " : ""
    logging && @info "Beginning classification of all $(primstring)$k-rational polygons with one interior lattice point."

    Ps = classify_maximal_polygons_genus_one(k; primitive, logging)
    return subpolygons(Ps; primitive, out_path, logging)
end

