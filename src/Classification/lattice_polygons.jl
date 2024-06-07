
@doc raw"""
    classify_maximal_lattice_polygons_with_collinear_interior_points(g :: Int, T :: Type{<:Integer} = Int)

Return all maximal lattice polygons with `g` collinear interior lattice
points.

"""
function classify_maximal_lattice_polygons_with_collinear_interior_points(g :: Int, T :: Type{<:Integer} = Int)
    Ps = RationalPolygon{T}[]
    g <= 1 && return Ps
    push!(Ps,RationalPolygon(LatticePoint{T}[(0,-1),(2g+2,-1),(0,1)], 1))
    for i = 1 : g+1
        push!(Ps,RationalPolygon(LatticePoint{T}[(0,-1),(2g+2-i,-1),(i,1),(0,1)], 1))
    end
    return Ps
end

@doc raw"""
    classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(g :: Int, T :: Type{<:Integer} = Int)

Return all maximal lattice polygons with `g` interior lattice points, where the
convex hull of these points is a two-dimensional lattice polygon without
interior lattice points.

"""
function classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(g :: Int, T :: Type{<:Integer} = Int)

    g <= 2 && return RationalPolygon{T}[]
    g == 3 && return [RationalPolygon(LatticePoint{T}[(0,0),(4,0),(0,4)], 1)]

    Ps = RationalPolygon{T}[]

    if g % 3 == 1
        push!(Ps, RationalPolygon(LatticePoint{T}[(-1,-1),(g+1,-1),(-1,2)], 1))
    end
    imin = g % 3 == 1 ? ((g+1) ÷ 3) + 1 : ((g+1) ÷ 3)
    imax = g ÷ 2
    for i = imin : imax
        push!(Ps, RationalPolygon(LatticePoint{T}[(-1,-1),(2g-3i,-1),(3i-g,2),(-1,2)], 1))
    end
    if g == 6
        push!(Ps, RationalPolygon(LatticePoint{T}[(0,0),(5,0),(0,5)], 1))
    end
    return Ps
end

function _is_internal_quick(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    for i = 1 : N
        u, v, w = scaled_vertex(P, i-1), scaled_vertex(P,i), scaled_vertex(P,i+1)
        p1, p2 = u - v, w - v
        q1, q2 = primitivize(p1), primitivize(p2)
        d, k, _ = cls_cone_normal_form(SA[q1[1] q2[1] ; q1[2] q2[2]])
        if (d,k) != (1,0) && k != 1
            return false
        end
    end
    return true
end

abstract type LatticePolygonsByGenusStorage{T <: Integer} end

mutable struct InMemoryLatticePolygonsByGenusStorage{T <: Integer} <: LatticePolygonsByGenusStorage{T}
    maximal_genus :: Int
    maximal_polygons :: Vector{Vector{RationalPolygon{T}}}
    all_polygons :: Vector{Vector{RationalPolygon{T}}}
    last_completed_genus :: Int
    total_count :: Int

    function InMemoryLatticePolygonsByGenusStorage{T}(maximal_genus :: Int) where {T <: Integer}
        all_polygons = Vector{Vector{RationalPolygon{T}}}(undef, maximal_genus)
        maximal_polygons = Vector{RationalPolygon{T}}[]
        push!(maximal_polygons, classify_maximal_polygons_genus_one(1))
        for i = 2 : maximal_genus
            push!(maximal_polygons, RationalPolygon{T}[])
        end

        return new{T}(maximal_genus, maximal_polygons, all_polygons, 0, 0)
    end

end

function classify_next_genus(st :: InMemoryLatticePolygonsByGenusStorage{T}) where {T <: Integer}

    i = st.last_completed_genus + 1

    Ps_max = st.maximal_polygons[i]
    append!(Ps_max, classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(i,T))
    append!(Ps_max, classify_maximal_lattice_polygons_with_collinear_interior_points(i,T))

    st.all_polygons[i] = subpolygons(Ps_max; normal_form = :affine)
    st.total_count += length(st.all_polygons[i])

    out_array = Vector{Tuple{Int,RationalPolygon{T}}}[]
    for i = 1 : Threads.nthreads()
        push!(out_array, Tuple{Int,RationalPolygon{T}}[])
    end

    Threads.@threads for P ∈ st.all_polygons[i]
        _is_internal_quick(P) || continue
        n = number_of_lattice_points(P)
        n ≤ st.maximal_genus || continue
        Q = move_out_edges(P)
        rationality(Q) == 1 || continue
        push!(out_array[Threads.threadid()], (n,Q))
    end
    moved_out_polygons = vcat(out_array...)

    for (n, P) ∈ moved_out_polygons
        push!(st.maximal_polygons[n], P)
    end

    st.last_completed_genus = i

    return length(st.all_polygons[i]), length(st.maximal_polygons[i]), st.total_count

end


mutable struct OnDiskLatticePolygonsByGenusStorage{T <: Integer} <: LatticePolygonsByGenusStorage{T}
    maximal_genus :: Int
    directory :: String
    last_completed_genus :: Int
    total_count :: Int

    function OnDiskLatticePolygonsByGenusStorage{T}(maximal_genus :: Int, directory :: String) where {T <: Integer}

        isdir(directory) || error("$directory is not a directory")
        isempty(readdir(directory)) || error("$directory is dirty. Please provide an empty directory")
        mkdir(joinpath(directory, "maximal"))
        mkdir(joinpath(directory, "all"))

        write_rational_polygons(classify_maximal_polygons_genus_one(one(T)), joinpath(directory, "maximal/i1.txt"))
        mkdir(joinpath(directory, "all/i1"))

        for i = 2 : maximal_genus
            touch(joinpath(directory, "maximal/i$i.txt"))
            mkdir(joinpath(directory, "all/i$i"))
        end

        new{T}(maximal_genus, directory, 0, 0)

    end

end

function classify_next_genus(st :: OnDiskLatticePolygonsByGenusStorage{T}) where {T <: Integer}

    i = st.last_completed_genus + 1
    max_file_path = joinpath(st.directory, "maximal/i$i.txt")
    all_path = joinpath(st.directory, "all/i$i")

    Ps_max = parse_rational_polygons(one(T), max_file_path) 
    append!(Ps_max, classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(i,T))
    append!(Ps_max, classify_maximal_lattice_polygons_with_collinear_interior_points(i,T))

    subpolygons(Ps_max; normal_form = :affine, out_path = all_path)

    out_array = Vector{Tuple{Int,RationalPolygon{T}}}[]
    for i = 1 : Threads.nthreads()
        push!(out_array, Tuple{Int,RationalPolygon{T}}[])
    end

    files = readdir(all_path)
    filter!(f -> startswith(f, "vol_"), files)
    all_count = 0
    for f ∈ files
        Ps = parse_rational_polygons(one(T), joinpath(all_path, f))
        all_count += length(Ps)

        for P ∈ Ps
            _is_internal_quick(P) || continue
            n = number_of_lattice_points(P)
            n ≤ st.maximal_genus || continue
            Q = move_out_edges(P)
            rationality(Q) == 1 || continue
            push!(out_array[Threads.threadid()], (n,Q))
        end
    end
    st.total_count += all_count

    moved_out_polygons = vcat(out_array...)

    for (n, P) ∈ moved_out_polygons
        filepath = joinpath(st.directory, "maximal/i$n.txt")
        write_rational_polygons([P], filepath, "a")
    end

    write_rational_polygons(Ps_max, max_file_path, "w")

    st.last_completed_genus = i

    return all_count, length(Ps_max), st.total_count

end


lattice_polygons_by_genus_storage(maximal_genus :: Int, directory :: Union{Missing,String} = missing, T :: Type{<:Integer} = Int) =
ismissing(directory) ? InMemoryLatticePolygonsByGenusStorage{T}(maximal_genus) : OnDiskLatticePolygonsByGenusStorage{T}(maximal_genus, directory)


function classify_lattice_polygons_by_genus(
        st :: LatticePolygonsByGenusStorage{T}; 
        logging :: Bool = false) where {T <: Integer}

    for i = st.last_completed_genus + 1 : st.maximal_genus

        all_count, max_count, total_count = classify_next_genus(st)

        logging && @info "Genus $i. Number of (maximal) polygons: $all_count ($max_count). Total: $total_count."

    end

    if st isa InMemoryLatticePolygonsByGenusStorage
        return (st.all_polygons, st.maximal_polygons)
    else
        return st
    end

end

@doc raw"""
    classify_lattice_polygons_by_genus(g :: Int, T :: Type{<:Integer} = Int; out_path :: Union{Missing,String} = missing, logging :: Bool = false)

Classify all lattice polygons with up to `g` interior lattice points. If
`out_path` is specified, the polygons will be saved to that location, with one
text file for each number `1 ≤ i ≤ g` of interior lattice points. If `out_path`
is not specified, the resulting polygons will be kept in memory and returned as
a tuple `(Ps, Ps_max)`, where `Ps` and `Ps_max` are vectors of length `g`
containing all polygons and just the maximal polygons respectively.

"""
classify_lattice_polygons_by_genus(
        g :: Int, T :: Type{<:Integer} = Int; 
        out_path :: Union{Missing,String} = missing,
        logging :: Bool = false) =
classify_lattice_polygons_by_genus(lattice_polygons_by_genus_storage(g, out_path, T); logging)


@doc raw"""
    filter_fano_polygons(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}   
Given a list of rational polygons `Ps`, return the list of all fano polygons
that are affine equivalent to a polygons from `Ps`.

"""
function filter_fano_polygons(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    res = RationalPolygon{T}[]
    for P ∈ Ps, p ∈ interior_lattice_points(P)
        Q = P - p
        is_fano(Q) || continue
        all(Q2 -> !are_unimodular_equivalent(Q,Q2), res) || continue
        push!(res,Q)
    end
    return res
end
