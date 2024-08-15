
@doc raw"""
    classify_maximal_polygons_m1p1(k :: T, i :: Int)

Return all maximal `k`-rational polygons with `i` interior lattice points
that can be realized in $\mathbb{R} \times [-1,1]$.

"""
function classify_maximal_polygons_m1p1(k :: T, i :: Int) where {T <: Integer}
    A = RationalPolygon(SA[0 k 0 ; 0 k k], k)
    B = RationalPolygon(SA[0 k*(i+1) k*(i+2) 0 ; 0 0 k k], k)

    vs = filter(v -> v[2] > 0, k_rational_points(A,k))
    ws = filter(w -> w[2] > 0, k_rational_points(B,k))

    Pss = Set{RationalPolygon{T}}[]
    for k = 1 : Threads.nthreads()
        push!(Pss, Set{RationalPolygon{T}}())
    end

    Threads.@threads for v ∈ vs
        tid = Threads.threadid()
        for w ∈ ws
            H1 = affine_halfplane(v,RationalPoint{T}(0,0))
            H2 = affine_halfplane(RationalPoint{T}(i+1,0),w)
            w ∈ H1 && v ∈ H2 || continue

            H_upper = affine_halfplane(RationalPoint{T}(0,-1),-T(1))
            H_lower = affine_halfplane(RationalPoint{T}(0,1),-T(1))

            P = k_rational_hull(intersect_halfplanes([H1,H2,H_upper,H_lower]), k)
            number_of_interior_lattice_points(P) == i || continue
            P ∉ Pss[tid] || continue
            is_maximal(P) || continue

            push!(Pss[tid], affine_normal_form(P))
        end
    end

    return union(Pss...)
end
