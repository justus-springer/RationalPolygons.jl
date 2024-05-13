
@doc raw"""
    remove_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

Remove the `i`-th vertex from `P`, i.e. return the convex hull of all
`k`-rational points of `P` except the `i`-th vertex.

"""
function remove_vertex(P :: RationalPolygon{T,N}, i :: Int) where {N,T <: Integer}
    u, v, w = lattice_vertex(P, i-1), lattice_vertex(P,i), lattice_vertex(P,i+1)
    p1, p2 = u - v, w - v
    q1, q2 = primitivize(p1), primitivize(p2)
    hb = [p + v for p ∈ hilbert_basis(q1,q2)]

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

    return Q

end



@doc raw"""
    subpolygons(starting_polygons :: Vector{<:RationalPolygon{T}}; out_path :: Union{Missing, String} = missing) where {T <: Integer}

Given some rational polygons `starting_polygons` with shared rationality `k`
and number of interior lattice points `n`, compute all subpolygons sharing
the same rationality and number of interior lattice points.

"""
function subpolygons(starting_polygons :: Vector{<:RationalPolygon{T}}; out_path :: Union{Missing, String} = missing, logging = false) where {T <: Integer}

    logging && @info "Starting to compute subpolygons..."

    k = rationality(first(starting_polygons))
    n = number_of_interior_lattice_points(first(starting_polygons))
    all(P -> rationality(P) == k, starting_polygons) || error("all polygons must have the same rationality")
    all(P -> number_of_interior_lattice_points(P) == n, starting_polygons) || error("all polygons must have the same number of interior lattice points")

    st = subpolygon_storage(starting_polygons, out_path)
    Ps, current_area = next_polygons!(st)
    total_count = length(Ps)

    while current_area > 3

        logging && @info "Volume: $current_area. Number of polygons: $(length(Ps)). Total: $total_count"

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
                    Q = normal_form(remove_vertex(P,j))
                    number_of_vertices(Q) > 2 || continue
                    number_of_interior_lattice_points(Q) == n || continue
                    push!(out_array[Threads.threadid()], Q)
                end
            end
        end
        save!(st, vcat(out_array...))
        
        Ps, current_area = next_polygons!(st)
        total_count += length(Ps)
        
    end

    logging && @info "Volume: $current_area. Number of polygons: $(length(Ps)). Total: $total_count"

    if ismissing(out_path)
        return Ps = vcat(values(st.polygons_dict)...)
    else
        return nothing
    end
                    
end

subpolygons(P :: RationalPolygon{T}) where {T <: Integer} =
subpolygons([P])
