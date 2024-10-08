
@doc raw"""
    classify_maximal_polygons_genus_one_m1p1(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{R} \times [-1,1]$. If `primitive = true`
is passed, then only primitive polygons (i.e. ldp polygons) are returned.

"""
function classify_maximal_polygons_genus_one_m1p1(k :: T; primitive :: Bool = false) where {T <: Integer}

    A = RationalPolygon(SA[-k 0 -k ; 0 k k], k)
    B = RationalPolygon(SA[-k k 2k -k ; 0 0 k k], k)

    vs = filter(v -> v[2] > 0, k_rational_points(A,k))
    primitive && filter!(v -> is_primitive(k*v), vs)
    ws = filter(w -> w[2] > 0, k_rational_points(B,k))
    primitive && filter!(w -> is_primitive(k*w), ws)

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v ∈ vs
        tid = Threads.threadid()
        for w ∈ ws
            H1 = affine_halfplane(v,RationalPoint{T}(-1,0))
            H2 = affine_halfplane(RationalPoint{T}(1,0),w)
            w ∈ H1 && v ∈ H2 || continue

            H_upper = affine_halfplane(RationalPoint{T}(0,-1),-T(1))
            H_lower = affine_halfplane(RationalPoint{T}(0,1),-T(1))

            P = k_rational_hull(intersect_halfplanes([H1,H2,H_upper,H_lower]), k; primitive)
            interior_lattice_points(P) == [zero(LatticePoint{T})] || continue
            P ∉ Pss[tid] || continue

            if !primitive
                is_maximal(P) || continue
            end

            push!(Pss[tid], unimodular_normal_form(P))
        end
    end

    return union!(Pss...)
end




@doc raw"""
    classify_maximal_polygons_genus_one_m1p2(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{R} \times [-1,2]$. If `primitive = true`
is passed, then only primitive polygons (i.e. ldp polygons) are returned.

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

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v1 ∈ vs
        tid = Threads.threadid()
        for v2 ∈ vs
            Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
            v1 ∈ Ha2 && v2 ∈ Ha1 || continue

            x1 = intersection_point(line(Ha1), horizontal_line(0))[1]
            x2 = intersection_point(line(Ha2), horizontal_line(0))[1]
            ws1_or_nothing = x1 > -1 ? [zero(RationalPoint{T})] : ws
            ws2_or_nothing = x2 < 1 ? [zero(RationalPoint{T})] : ws

            for (w1,w2) ∈ Iterators.product(ws1_or_nothing, ws2_or_nothing)

                w1_is_nothing = w1 == zero(RationalPoint{T})
                w2_is_nothing = w2 == zero(RationalPoint{T})

                w1_is_nothing || (w1 ∈ Ha1 && w1 ∈ Ha2 || continue)
                w2_is_nothing || (w2 ∈ Ha1 && w2 ∈ Ha2 || continue)

                H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
                H_lower = affine_halfplane(RationalPoint{T}(0,-1),RationalPoint{T}(1,-1))
                halfplanes = [Ha1, Ha2, H_upper, H_lower]
                
                if !w1_is_nothing
                    Hb1 = affine_halfplane(b1,w1)
                    push!(halfplanes, Hb1)
                    v1 ∈ Hb1 && v2 ∈ Hb1 || continue
                    w2_is_nothing || w2 ∈ Hb1 || continue
                end

                if !w2_is_nothing
                    Hb2 = affine_halfplane(w2,b2)
                    push!(halfplanes, Hb2)
                    v1 ∈ Hb2 && v2 ∈ Hb2 || continue
                    w1_is_nothing || w1 ∈ Hb2 || continue
                end

                all(contains_origin_in_interior, halfplanes) || continue

                P = k_rational_hull(intersect_halfplanes(halfplanes), k; primitive)

                !primitive && (is_maximal(P) || continue)

                push!(Pss[tid], unimodular_normal_form(P))
            end
        end
    end

    return union(Pss...)

end

@doc raw"""
    classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point that can be realized in $\mathbb{R} \times [-2,2]$ that have non-empty
intersection with the `q`-th classification box, where `1 ≤ q ≤ 3`. If
`primitive = true` is passed, then only primitive polygons (i.e. ldp polygons)
are returned.

"""
function classify_maximal_polygons_genus_one_m2p2(k :: T, q :: Int; primitive :: Bool = false) where {T <: Integer}
    

    A = RationalPolygon(SA[-2k -k 0 0 ; 2k k k 2k], k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(1,0)

    if q == 1
        B1 = RationalPolygon(SA[-2k -k -k ; -k -k 0], k)
        B2 = RationalPolygon(SA[-2k k k ; -k -k 0], k)
        C = RationalPolygon(SA[-3k -2k -k -2k ; -2k -2k -k -k], k)
        c1,c2 = RationalPoint{T}(-2,-1), RationalPoint{T}(-1,-1)
    elseif q == 2
        B1 = RationalPolygon(SA[-2k 0 -k ; -k -k 0], k)
        B2 = RationalPolygon(SA[-k 2k k ; -k -k 0], k)
        C = RationalPolygon(SA[-2k 0 0 -k ; -2k -2k -k -k], k)
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

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v1 ∈ vs
        tid = Threads.threadid()
        for v2 ∈ vs
            Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
            v1 ∈ Ha2 && v2 ∈ Ha1 || continue

            x1 = intersection_point(line(Ha1), horizontal_line(0))[1]
            x2 = intersection_point(line(Ha2), horizontal_line(0))[1]
            ws1_or_nothing = x1 > b1[1] ? [zero(RationalPoint{T})] : ws1
            ws2_or_nothing = x2 < b2[1] ? [zero(RationalPoint{T})] : ws2

            y1 = intersection_point(line(Ha1), horizontal_line(-1))[1]
            y2 = intersection_point(line(Ha2), horizontal_line(-1))[1]
            y2 - y1 ≥ 3 || continue
            
            for (w1,w2) ∈ Iterators.product(ws1_or_nothing, ws2_or_nothing)
                
                w1_is_nothing = w1 == zero(RationalPoint{T})
                w2_is_nothing = w2 == zero(RationalPoint{T})

                w1_is_nothing || (w1 ∈ Ha1 && w1 ∈ Ha2 || continue)
                w2_is_nothing || (w2 ∈ Ha1 && w2 ∈ Ha2 || continue)

                H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
                H_lower = affine_halfplane(RationalPoint{T}(0,-2),RationalPoint{T}(1,-2))
                halfplanes = [Ha1, Ha2, H_upper, H_lower]
                
                if !w1_is_nothing
                    Hb1 = affine_halfplane(b1,w1)
                    push!(halfplanes, Hb1)
                    v1 ∈ Hb1 && v2 ∈ Hb1 || continue
                    w2_is_nothing || w2 ∈ Hb1 || continue
                    y1_updated = max(y1, intersection_point(line(Hb1), horizontal_line(-1))[1])
                else
                    y1_updated = y1
                end

                if !w2_is_nothing
                    Hb2 = affine_halfplane(w2,b2)
                    push!(halfplanes, Hb2)
                    v1 ∈ Hb2 && v2 ∈ Hb2 || continue
                    w1_is_nothing || w1 ∈ Hb2 || continue
                    y2_updated = min(y2, intersection_point(line(Hb2), horizontal_line(-1))[1])
                else
                    y2_updated = y2
                end

                us1_or_nothing = y1_updated > c1[1] ? [zero(RationalPoint{T})] : us
                us2_or_nothing = y2_updated < c2[1] ? [zero(RationalPoint{T})] : us

                for (u1,u2) ∈ Iterators.product(us1_or_nothing, us2_or_nothing)

                    u1_is_nothing = u1 == zero(RationalPoint{T})
                    u2_is_nothing = u2 == zero(RationalPoint{T})

                    u1_is_nothing || (u1 ∈ Ha1 && u1 ∈ Ha2 || continue)
                    u2_is_nothing || (u2 ∈ Ha1 && u2 ∈ Ha2 || continue)

                    updated_halfplanes = copy(halfplanes)

                    if !u1_is_nothing
                        Hc1 = affine_halfplane(c1,u1)
                        push!(updated_halfplanes, Hc1)
                        v1 ∈ Hc1 && v2 ∈ Hc1 || continue
                        w1_is_nothing || (w1 ∈ Hc1 && u1 ∈ Hb1) || continue
                        w2_is_nothing || (w2 ∈ Hc1 && u1 ∈ Hb2) || continue
                        u2_is_nothing || u2 ∈ Hc1 || continue
                    end

                    if !u2_is_nothing
                        Hc2 = affine_halfplane(u2,c2)
                        push!(updated_halfplanes, Hc2)
                        v1 ∈ Hc2 && v2 ∈ Hc2 || continue
                        w1_is_nothing || (w1 ∈ Hc2 && u2 ∈ Hb1) || continue
                        w2_is_nothing || (w2 ∈ Hc2 && u2 ∈ Hb2) || continue
                        u1_is_nothing || u1 ∈ Hc2 || continue
                    end

                    all(contains_origin_in_interior, updated_halfplanes) || continue

                    P = k_rational_hull(intersect_halfplanes(updated_halfplanes), k; primitive)

                    !primitive && (is_maximal(P) || continue)

                    push!(Pss[tid], unimodular_normal_form(P))
                end
            end
        end
    end

    return union!(Pss...)

end


@doc raw"""
    classify_maximal_polygons_genus_one(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with exactly one interior lattice
point. If `primitive = true` is passed, then only primitive polygons (i.e. ldp
polygons) are returned.

# Example

Compute the numbers of polygons for ``k ≤ 3``.

```jldoctest
julia> length.(classify_maximal_polygons_genus_one.(1:3))
4-element Vector{Int64}:
   3
  10
  39
```

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

- `primitive :: Bool`. If set to true, only primitive polygons (i.e. ldp polygons) are returned.
- `logging :: Bool`. Controls whether to display logging messages showing the
current progress.

# Example

Reproduce the classifcation of all 5145 half-integral polygons with exactly one
interior lattice point. It first computes all maximal polygons with
`classify_maximal_polygons_genus_one` and then generates all their subpolygons.

```jldoctest
julia> classify_polygons_genus_one(2; logging=true);
[ Info: Found 9 maximal polygons in QQ x [-1,1]. New: 9, total: 9
[ Info: Found 2 maximal polygons in QQ x [-1,2]. New: 1, total: 10
[ Info: Found 0 maximal polygons in QQ x [-2,2], box 1. New: 0, total: 10
[ Info: Found 3 maximal polygons in QQ x [-2,2], box 2. New: 0, total: 10
[ Info: Found 3 maximal polygons in QQ x [-2,2], box 3. New: 0, total: 10
[ Info: [a = 36]. Polygons to peel: 2.
[ Info: [a = 36]. Peeling complete. New polygons: 2. Running total: 12
[ Info: [a = 35]. Polygons to peel: 2.
[ Info: [a = 35]. Peeling complete. New polygons: 5. Running total: 17
[ Info: [a = 34]. Polygons to peel: 5.
[ Info: [a = 34]. Peeling complete. New polygons: 10. Running total: 27
[ Info: [a = 33]. Polygons to peel: 9.
[ Info: [a = 33]. Peeling complete. New polygons: 19. Running total: 46
[ Info: [a = 32]. Polygons to peel: 23.
[ Info: [a = 32]. Peeling complete. New polygons: 45. Running total: 91
[ Info: [a = 31]. Polygons to peel: 31.
[ Info: [a = 31]. Peeling complete. New polygons: 69. Running total: 160
[ Info: [a = 30]. Polygons to peel: 60.
[ Info: [a = 30]. Peeling complete. New polygons: 115. Running total: 275
[ Info: [a = 29]. Polygons to peel: 84.
[ Info: [a = 29]. Peeling complete. New polygons: 170. Running total: 445
[ Info: [a = 28]. Polygons to peel: 137.
[ Info: [a = 28]. Peeling complete. New polygons: 239. Running total: 684
[ Info: [a = 27]. Polygons to peel: 171.
[ Info: [a = 27]. Peeling complete. New polygons: 285. Running total: 969
[ Info: [a = 26]. Polygons to peel: 240.
[ Info: [a = 26]. Peeling complete. New polygons: 364. Running total: 1333
[ Info: [a = 25]. Polygons to peel: 286.
[ Info: [a = 25]. Peeling complete. New polygons: 440. Running total: 1773
[ Info: [a = 24]. Polygons to peel: 356.
[ Info: [a = 24]. Peeling complete. New polygons: 466. Running total: 2239
[ Info: [a = 23]. Polygons to peel: 351.
[ Info: [a = 23]. Peeling complete. New polygons: 433. Running total: 2672
[ Info: [a = 22]. Polygons to peel: 396.
[ Info: [a = 22]. Peeling complete. New polygons: 439. Running total: 3111
[ Info: [a = 21]. Polygons to peel: 391.
[ Info: [a = 21]. Peeling complete. New polygons: 394. Running total: 3505
[ Info: [a = 20]. Polygons to peel: 390.
[ Info: [a = 20]. Peeling complete. New polygons: 338. Running total: 3843
[ Info: [a = 19]. Polygons to peel: 349.
[ Info: [a = 19]. Peeling complete. New polygons: 269. Running total: 4112
[ Info: [a = 18]. Polygons to peel: 357.
[ Info: [a = 18]. Peeling complete. New polygons: 249. Running total: 4361
[ Info: [a = 17]. Polygons to peel: 301.
[ Info: [a = 17]. Peeling complete. New polygons: 215. Running total: 4576
[ Info: [a = 16]. Polygons to peel: 292.
[ Info: [a = 16]. Peeling complete. New polygons: 178. Running total: 4754
[ Info: [a = 15]. Polygons to peel: 233.
[ Info: [a = 15]. Peeling complete. New polygons: 121. Running total: 4875
[ Info: [a = 14]. Polygons to peel: 191.
[ Info: [a = 14]. Peeling complete. New polygons: 91. Running total: 4966
[ Info: [a = 13]. Polygons to peel: 140.
[ Info: [a = 13]. Peeling complete. New polygons: 71. Running total: 5037
[ Info: [a = 12]. Polygons to peel: 121.
[ Info: [a = 12]. Peeling complete. New polygons: 40. Running total: 5077
[ Info: [a = 11]. Polygons to peel: 67.
[ Info: [a = 11]. Peeling complete. New polygons: 21. Running total: 5098
[ Info: [a = 10]. Polygons to peel: 56.
[ Info: [a = 10]. Peeling complete. New polygons: 14. Running total: 5112
[ Info: [a = 9]. Polygons to peel: 38.
[ Info: [a = 9]. Peeling complete. New polygons: 13. Running total: 5125
[ Info: [a = 8]. Polygons to peel: 31.
[ Info: [a = 8]. Peeling complete. New polygons: 13. Running total: 5138
[ Info: [a = 7]. Polygons to peel: 14.
[ Info: [a = 7]. Peeling complete. New polygons: 3. Running total: 5141
[ Info: [a = 6]. Polygons to peel: 13.
[ Info: [a = 6]. Peeling complete. New polygons: 3. Running total: 5144
[ Info: [a = 5]. Polygons to peel: 4.
[ Info: [a = 5]. Peeling complete. New polygons: 1. Running total: 5145
[ Info: [a = 4]. Polygons to peel: 3.
[ Info: [a = 4]. Peeling complete. New polygons: 0. Running total: 5145
[ Info: [a = 3]. Polygons to peel: 1.
[ Info: [a = 3]. Peeling complete. New polygons: 0. Running total: 5145
[ Info: [a = 2]. Polygons to peel: 0.
[ Info: [a = 2]. Peeling complete. New polygons: 0. Running total: 5145
[ Info: [a = 1]. Polygons to peel: 0.
[ Info: [a = 1]. Peeling complete. New polygons: 0. Running total: 5145
```

"""
classify_polygons_genus_one(k :: T; primitive :: Bool = false, logging = false) where {T <: Integer} =
subpolygons(classify_maximal_polygons_genus_one(k; primitive, logging);
    use_affine_normal_form = false,
    only_equal_number_of_interior_lattice_points = true,
    primitive,
    logging)

