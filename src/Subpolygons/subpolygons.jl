
@doc raw"""
    remove_vertex(P :: RationalPolygon{T}, i :: Int) where {T <: Integer}

Remove the `i`-th vertex from `P`, i.e. return the convex hull of all
`k`-rational points of `P` except the `i`-th vertex.

"""
function remove_vertex(P :: RationalPolygon{T,N}, i :: Int; primitive :: Bool = false) where {N,T <: Integer}

    k = rationality(P)
    keeps_genus :: Bool = false

    if primitive
        v0 = vertex(P,i)
        ps = filter(p -> is_primitive(k*p), k_rational_points(P, k))
        filter!(p -> p != v0, ps)
        Q = convex_hull(ps, k)
        keeps_genus = number_of_interior_lattice_points(P) == number_of_interior_lattice_points(Q)
        return (Q, keeps_genus)
    else
        # the shortcut via hilbert bases only works in the non-primitive
        # case at the moment
        u, v, w = scaled_vertex(P, i-1), scaled_vertex(P,i), scaled_vertex(P,i+1)
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
            push!(vs, scaled_vertex(P,j))
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
    subpolygons(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    subpolygons(P :: RationalPolygon{T}) where {T <: Integer}
    subpolygons(st :: SubpolygonStorage{T}) where {T <: Integer}

Given a list of rational polygons with shared rationality and number of
interior lattice points, compute all subpolygons with the same number of
)nterior lattice points. The following keyword arguments are supported:

- `primitive :: Bool`: If set to true, only subpolygons with primitive vertices
are returned.

- `normal_form :: Symbol`: Can be either `:unimodular` or `:affine`. Used to
control which normal form should be used when comparing subpolygons for
equivalence.

- `out_path :: Union{Missing, String}`: If out_path is `missing`, then all
polygons are kept in memory. By specifying `out_path` to be a path to an empty
directory, the storage of the resulting polygons is delegated to the disk. In
this case, `subpolygons` will save the polygons to text files according to
their normalized area. Additionally, a file `last_volume` will be created and
constantly updated to hold the last fully completed area. This is useful for
resuming a lengthy computation at a later point, see also
`OnDiskSubpolygonStorage` for details on how to resume an unfinished
computations of subpolygons.

- `logging :: Bool`: Controls whether to show log messages showing the
progress.

"""
function subpolygons(Ps :: Vector{<:RationalPolygon{T}}; 
        primitive :: Bool = false, 
        normal_form :: Symbol = :unimodular,
        out_path :: Union{Missing, String} = missing, 
        logging :: Bool = false) where {T <: Integer}

    st = subpolygon_storage(Ps, out_path)
    return subpolygons(st; primitive, normal_form, logging)
                    
end

subpolygons(P :: RationalPolygon{T},
        primitive :: Bool = false, 
        normal_form :: Symbol = :unimodular,
        out_path :: Union{Missing, String} = missing, 
        logging :: Bool = false) where {T <: Integer} =
subpolygons([P]; primitive, normal_form, out_path, logging)
