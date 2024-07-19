
@doc raw"""
    classify_maximal_lattice_free_polygons_m1p1(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with no interior lattice points that
can be realized in $\mathbb{Q} \times [-1,1]$.

"""
function classify_maximal_lattice_free_polygons_m1p1(k :: T) where {T <: Integer}

    A = convex_hull(RationalPoint{T}.([(-1,0),(0,1),(-1,1)]), k)
    B = convex_hull(RationalPoint{T}.([(-1,0),(0,0),(1,1),(-1,1)]), k)

    vs = [v for v ∈ k_rational_points(A,k) if v[2] > 0]
    ws = [w for w ∈ k_rational_points(B,k) if w[2] > 0]

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v ∈ vs
        tid = Threads.threadid()
        for w ∈ ws
            H1 = affine_halfplane(v,RationalPoint{T}(-1,0))
            H2 = affine_halfplane(RationalPoint{T}(0,0),w)
            w ∈ H1 && v ∈ H2 || continue

            H_upper = affine_halfplane(RationalPoint{T}(0,-1),-T(1))
            H_lower = affine_halfplane(RationalPoint{T}(0,1),-T(1))

            P = k_rational_hull(intersect_halfplanes([H1,H2,H_upper,H_lower]), k)
            number_of_interior_lattice_points(P) == 0 || continue
            P ∉ Pss[tid] || continue
            is_maximal(P) || continue

            push!(Pss[tid], affine_normal_form(P))
        end
    end

    return union!(Pss...)
    
end


function classify_maximal_lattice_free_polygons_m1p2_squares(k :: T) where {T <: Integer}

    A = RationalPolygon(SMatrix{2,4,T}(0,k,k,k,k,2k,0,2k), k)
    B = RationalPolygon(SMatrix{2,4,T}(0,0,0,-k,k,-k,k,0), k)

    a1,a2 = RationalPoint{T}(0,1), RationalPoint{T}(1,1)
    b1,b2 = RationalPoint{T}(0,0), RationalPoint{T}(1,0)
    
    vs = filter(v -> v[2] > 1, k_rational_points(A, k))
    ws = filter(w -> w[2] < 0, k_rational_points(B, k))

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v1 ∈ vs
        tid = Threads.threadid()
        is_primitive(k*(v1 - a1)) || continue
        for v2 ∈ vs
            is_primitive(k*(v2 - a2)) || continue

            Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
            v1 ∈ Ha2 && v2 ∈ Ha1 || continue

            dv1, dv2 = direction_vector(Ha1), direction_vector(Ha2)
            dv1[2]*dv2[1] ≤ -dv2[2]*dv1[1] || continue

            for w1 ∈ ws, w2 ∈ ws
                is_primitive(k*(w1 - b1)) || continue
                is_primitive(k*(w2 - b2)) || continue

                Hb1, Hb2 = affine_halfplane(b1,w1), affine_halfplane(w2,b2)
                w1 ∈ Hb2 && w2 ∈ Hb1 || continue

                H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
                H_lower = affine_halfplane(RationalPoint{T}(0,-1),RationalPoint{T}(1,-1))

                P = k_rational_hull(intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]), k)
                number_of_vertices(P) > 2 || continue
                number_of_interior_lattice_points(P) == 0 || continue
                P ∉ Pss[tid] || continue
                is_maximal(P) || continue

                push!(Pss[tid], affine_normal_form(P))
            end
        end
    end

    return union!(Pss...)
end


function classify_maximal_lattice_free_polygons_m1p2_trapezoids(k :: T) where {T <: Integer}

    A1 = RationalPolygon(SMatrix{2,3,T}(-k,2k,0,k,0,2k), k)
    A2 = RationalPolygon(SMatrix{2,4,T}(-k,2k,0,k,k,k,k,2k), k)
    B = RationalPolygon(SMatrix{2,4,T}(0,0,0,-k,2k,-k,k,0), k)

    a1,a2 = RationalPoint{T}(0,1), RationalPoint{T}(1,1)
    b = RationalPoint{T}(1,0)
    
    vs1 = filter(v -> v[2] > 1, k_rational_points(A1, k))
    vs2 = filter(v -> v[2] > 1, k_rational_points(A2, k))
    ws = filter(w -> w[2] < 0, k_rational_points(B, k))

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v1 ∈ vs1
        is_primitive(k*(v1 - a1)) || continue
        tid = Threads.threadid()
        for v2 ∈ vs2
            is_primitive(k*(v2 - a2)) || continue
            Ha1, Ha2 = affine_halfplane(v1,a1), affine_halfplane(a2,v2)
            v1 ∈ Ha2 && v2 ∈ Ha1 || continue

            for w ∈ ws
                is_primitive(k*(w - b)) || continue
                w ∈ Ha1 && w ∈ Ha2 || continue

                Hb = affine_halfplane(w,b)
                v1 ∈ Hb && v2 ∈ Hb || continue

                H_upper = affine_halfplane(RationalPoint{T}(1,2),RationalPoint{T}(0,2))
                H_lower = affine_halfplane(RationalPoint{T}(0,-1),RationalPoint{T}(1,-1))

                P = k_rational_hull(intersect_halfplanes([Ha1,Ha2,Hb,H_upper,H_lower]), k)
                number_of_vertices(P) > 2 || continue
                number_of_interior_lattice_points(P) == 0 || continue
                P ∉ Pss[tid] || continue
                is_maximal(P) || continue

                push!(Pss[tid], affine_normal_form(P))
            end
        end
    end

    return union!(Pss...)
end


@doc raw"""
    classify_maximal_lattice_free_polygons(k :: T ; logging = false) where {T <: Integer}

Return all `k`-rational polygons with no interior lattice points.

"""
function classify_maximal_lattice_free_polygons(k :: T ; logging = false) where {T <: Integer}
    Ps = RationalPolygon{T}[]
    count, total_count = 0, 0
    
    new_Ps = classify_maximal_lattice_free_polygons_m1p1(k)
    count = length(new_Ps)
    union!(Ps, new_Ps)
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count lattice-free maximal polygons in QQ x [-1,1]. New: $new_count, total: $total_count"

    new_Ps = classify_maximal_lattice_free_polygons_m1p2_squares(k)
    count = length(new_Ps)
    union!(Ps, new_Ps)
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count lattice-free maximal polygons in QQ x [-1,2], squares case. New: $new_count, total: $total_count"

    new_Ps = classify_maximal_lattice_free_polygons_m1p2_trapezoids(k)
    count = length(new_Ps)
    union!(Ps, new_Ps)
    new_count = length(Ps) - total_count
    total_count = length(Ps)

    logging && @info "Found $count lattice-free maximal polygons in QQ x [-1,2], trapezoids case. New: $new_count, total: $total_count"


    return collect(Ps)
end
