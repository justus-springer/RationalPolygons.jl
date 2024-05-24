
@doc raw"""
    remove_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

Remove the `i`-th vertex from `P`, i.e. return the convex hull of all
`k`-rational points of `P` except the `i`-th vertex.

"""
function remove_vertex(P :: RationalPolygon{T,N}, i :: Int; primitive :: Bool = false) where {N,T <: Integer}

    k = rationality(P)

    if primitive
        v = vertex(P,i)
        ps = k_rational_points(k, P; primitive)
        filter!(p -> p != v, ps)
        Q = convex_hull(ps, k)
        keeps_genus = number_of_interior_lattice_points(P) == number_of_interior_lattice_points(Q)
        return (Q, keeps_genus)
    else
        # the shortcut via hilbert bases only works in the non-primitive
        # case at the moment
        u, v, w = lattice_vertex(P, i-1), lattice_vertex(P,i), lattice_vertex(P,i+1)
        p1, p2 = u - v, w - v
        q1, q2 = primitivize(p1), primitivize(p2)
        hb = [p + v for p ∈ hilbert_basis(q1,q2)]
        
        keeps_genus = all(p -> gcd(p) % k != 0, hb[2 : end-1])

        vs = LatticePoint{T}[]
        !is_primitive(p1) && push!(vs, u)
        push!(vs,first(hb))
        for j = 2 : length(hb)-1
            if hb[j]-hb[j-1] != hb[j+1]-hb[j]
                push!(vs, hb[j])
            end
        end
        push!(vs,last(hb))
        !is_primitive(p2) && push!(vs, w)
        for j = i+2 : i+N-2
            push!(vs, lattice_vertex(P,j))
        end

        Q = RationalPolygon(vs, rationality(P))

        return (Q, keeps_genus)

    end

end

function subpolygons(st :: SubpolygonStorage{T};
        primitive :: Bool = false,
        normal_form :: Symbol = :unimodular,
        logging :: Bool = false) where {T <: Integer}
    
    normal_form ∈ [:unimodular, :affine] || error("`normal_form` must be either `:unimodular` or `:affine`")

    n = number_of_interior_lattice_points(st)

    while !is_finished(st)

        Ps, current_area = next_polygons(st)

        logging && @info "Volume: $current_area. Number of polygons: $(length(Ps)). Total: $(total_count(st))"

        N = length(Ps)
        num_blocks = 2 * Threads.nthreads()
        b = N ÷ num_blocks

        out_array = Vector{RationalPolygon{T}}[]
        for i = 1 : Threads.nthreads()
            push!(out_array, RationalPolygon{T}[])
        end

        Threads.@threads for k = 1 : num_blocks
            lower_bound = (k-1) * b + 1
            upper_bound = k == num_blocks ? N : k * b
            for i = lower_bound : upper_bound
                P = Ps[i]
                for j = 1 : number_of_vertices(P)
                    Q, keeps_genus = remove_vertex(P, j; primitive)
                    keeps_genus || continue
                    number_of_vertices(Q) > 2 || continue
                    if normal_form == :unimodular
                        Q = unimodular_normal_form(Q)
                    else
                        Q = affine_normal_form(Q)
                    end
                    push!(out_array[Threads.threadid()], Q)
                end
            end
        end
        save!(st, vcat(out_array...))
        mark_volume_completed(st, current_area)
        
    end

    logging && @info "Found a total of $(total_count(st)) subpolygons"

    return return_value(st)

end


@doc raw"""
    subpolygons(starting_polygons :: Vector{<:RationalPolygon{T}}; out_path :: Union{Missing, String} = missing) where {T <: Integer}

Given some rational polygons `starting_polygons` with shared rationality `k`
and number of interior lattice points `n`, compute all subpolygons sharing
the same rationality and number of interior lattice points.

"""
function subpolygons(starting_polygons :: Vector{<:RationalPolygon{T}}; 
        primitive :: Bool = false, 
        normal_form :: Symbol = :unimodular,
        out_path :: Union{Missing, String} = missing, 
        logging = false) where {T <: Integer}

    st = subpolygon_storage(starting_polygons, out_path)
    return subpolygons(st; primitive, normal_form, logging)
                    
end

subpolygons(P :: RationalPolygon{T}) where {T <: Integer} =
subpolygons([P])
