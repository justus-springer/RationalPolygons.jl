

@doc raw"""
    classify_maximal_lattice_free_polygons_m1p2_squares(k :: T) where {T <: Integer}

Return all `k`-maximal polygons with no interior lattice points that are
contained in ``A \cup \mathbb{R} \times [-1,2] \cup B`` where ``A`` is the
square with vertices ``(0,1),(1,1),(1,2),(0,2)`` and B is the square with vertices
``(0,0),(0,-1),(1,-1),(1,0)``.

"""
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


@doc raw"""
    classify_maximal_lattice_free_polygons_m1p2_trapezoids(k :: T) where {T <: Integer}

Return all `k`-maximal polygons with no interior lattice points that are
contained in ``A \cup \mathbb{R} \times [-1,2] \cup B`` where ``A`` is the
trapezoid with vertices ``(-1,2),(0,1),(1,1),(1,2)`` and ``B`` is the trapezoid
with vertices ``(0,0),(0,-1),(2,-1),(1,0)``, excluding the polygons from
`classify_maximal_lattice_free_polygons_m1p2_squares`.

"""
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
    classify_maximal_lattice_free_polygons_m1p2(k :: T) where {T <: Integer}

Return all `k`-maximal polygons with no interior lattice points contained in
``\mathbb{R} \times [-1,2]``. This is simply the union of
`classify_maximal_lattice_free_polygons_m1p2_squares` and
`classify_maximal_lattice_free_polygons_m1p2_trapezoids`.

"""
classify_maximal_lattice_free_polygons_m1p2(k :: T) where {T <: Integer} =
union(classify_maximal_lattice_free_polygons_m1p2_squares(k),
      classify_maximal_lattice_free_polygons_m1p2_trapezoids(k))


@doc raw"""
    classify_maximal_lattice_free_polygons(k :: T ; logging = false) where {T <: Integer}

Return all `k`-rational polygons with no interior lattice points.

# Example

Compute the numbers of polygons for `1 ≤ k ≤ 6`:

```jldoctest
julia> length.(classify_maximal_lattice_free_polygons.(1:6))
6-element Vector{Int64}:
   1
   4
  14
  39
 134
 299
```

"""
function classify_maximal_lattice_free_polygons(k :: T ; logging = false) where {T <: Integer}
    Ps = RationalPolygon{T}[]
    count, total_count = 0, 0
    
    new_Ps = classify_maximal_polygons_m1p1(k,0)
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
