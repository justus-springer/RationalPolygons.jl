
@doc raw"""
    mutable struct SubpolygonStorage{T<:Integer}   

A struct used by the function `subpolygons` to hold intermediate results.

"""
mutable struct SubpolygonStorage{T<:Integer}
    dir :: Union{String, Missing}
    dict :: Dict{T,Set{RationalPolygon{T}}}
    last_volume :: T

    function SubpolygonStorage{T}(starting_polygons :: Vector{RationalPolygon{T}}, dir :: Union{String, Missing} = missing) where {T <: Integer}
        !isempty(starting_polygons) || error("must provide a non-empty list of starting polygons")
        if !ismissing(dir)
            isdir(dir) || error("$dir is not a directory")
        end

        a = maximum(area.(starting_polygons))
        st = new{T}(dir, Dict{T,Vector{RationalPolygon{T}}}(), a + 1)
        save!(st, starting_polygons)
        return st
    end

end

function save!(st :: SubpolygonStorage{T}, Ps :: Vector{RationalPolygon{T}}) where {T <: Integer}
    for P ∈ Ps
        a = area(P)
        if !haskey(st.dict, a) 
            st.dict[a] = Set{RationalPolygon{T}}()
        end
        P ∉ st.dict[a] || continue
        push!(st.dict[a], P)
    end
end

function next_polygons!(st :: SubpolygonStorage{T}) where {T <: Integer}
    # Get the polygons with maximal area
    a = maximum(filter(b -> b < st.last_volume, keys(st.dict)))
    Ps = st.dict[a]

    # save them to a file and delete from our dictionary
    if !ismissing(st.dir)
        open(joinpath(st.dir, "vol_$a.txt"), "w") do f
            for P ∈ Ps
                V = vertex_matrix(P)
                println(f, [Vector(V[:,i]) for i = 1 : number_of_vertices(P)])
            end
        end
        delete!(st.dict, a)
    end

    st.last_volume = a

    return (collect(Ps), a)

end

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

    st = SubpolygonStorage{T}(starting_polygons, out_path)
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

    return st
                    
end

subpolygons(P :: RationalPolygon{T}) where {T <: Integer} =
subpolygons!([P])
