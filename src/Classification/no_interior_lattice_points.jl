
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

    return collect(union!(Pss...))
    
end


@doc raw"""
    classify_maximal_lattice_free_polygons_m1p2(k :: T) where {T <: Integer}

Return all maximal `k`-rational polygons with no interior lattice points that
can be realized in $\mathbb{Q} \times [-1,2]$.

"""
function classify_maximal_lattice_free_polygons_m1p2(k :: T) where {T <: Integer}

    A = convex_hull(RationalPoint{T}.([(0,1),(1,2),(-2,2),(-1,1)]), k)
    B = convex_hull(RationalPoint{T}.([(0,0),(1,-1),(-2,-1),(-1,0)]), k)

    a1,a2 = RationalPoint{T}(-1,1), RationalPoint{T}(0,1)
    b1,b2 = RationalPoint{T}(-1,0), RationalPoint{T}(0,0)
    
    vs = filter(v -> v[2] > 1, k_rational_points(A, k))
    ws = filter(w -> w[2] < 0, k_rational_points(B, k))

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v1 ∈ vs
        tid = Threads.threadid()
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

                P = k_rational_hull(intersect_halfplanes([Ha1,Ha2,Hb1,Hb2,H_upper,H_lower]), k)
                number_of_vertices(P) > 2 || continue
                number_of_interior_lattice_points(P) == 0 || continue
                P ∉ Pss[tid] || continue
                is_maximal(P) || continue

                push!(Pss[tid], affine_normal_form(P))
            end
        end
    end

    return collect(union!(Pss...))

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
