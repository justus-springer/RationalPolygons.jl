
@doc raw"""
    classify_maximal_lattice_polygons_with_collinear_interior_points(i :: Int, T :: Type{<:Integer} = Int)

Return all maximal lattice polygons with `i` collinear interior lattice
points.

"""
function classify_maximal_lattice_polygons_with_collinear_interior_points(i :: Int, T :: Type{<:Integer} = Int)
    Ps = RationalPolygon{T}[]
    i <= 1 && return Ps
    push!(Ps,RationalPolygon(LatticePoint{T}[(0,-1),(2i+2,-1),(0,1)], 1))
    for j = 1 : i+1
        push!(Ps,RationalPolygon(LatticePoint{T}[(0,-1),(2i+2-j,-1),(j,1),(0,1)], 1))
    end
    return Ps
end

@doc raw"""
    classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(i :: Int, T :: Type{<:Integer} = Int)

Return all maximal lattice polygons with `i` interior lattice points, where the
convex hull of these points is a two-dimensional lattice polygon without
interior lattice points.

"""
function classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(i :: Int, T :: Type{<:Integer} = Int)

    i <= 2 && return RationalPolygon{T}[]
    i == 3 && return [RationalPolygon(LatticePoint{T}[(0,0),(4,0),(0,4)], 1)]

    Ps = RationalPolygon{T}[]

    if i % 3 == 1
        push!(Ps, RationalPolygon(LatticePoint{T}[(-1,-1),(i+1,-1),(-1,2)], 1))
    end
    jmin = i % 3 == 1 ? ((i+1) ÷ 3) + 1 : ((i+1) ÷ 3)
    jmax = i ÷ 2
    for j = jmin : jmax
        push!(Ps, RationalPolygon(LatticePoint{T}[(-1,-1),(2i-3j,-1),(3j-i,2),(-1,2)], 1))
    end
    if i == 6
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

@doc raw"""
    abstract type CastryckStorage{T <: Integer} end   

Abstract supertype of `InMemoryCastryckStorage` and `HDFCastryckStorage`. Both
implement `classify_next_genus`, which performs a single step in Castryck's
classification of lattice polygons.

"""
abstract type CastryckStorage{T <: Integer} end


@doc raw"""
    mutable struct InMemoryCastryckStorage{T <: Integer} <: CastryckStorage{T}   

A struct holding classification results of Castryck's classification of lattice
polygons by number of interior lattice points (i.e. their genus). It has the following fields:

- `maximum_genus :: Int`: An upper bound for the maximum genus that the
   classification should run to. Defaults to `100`.
- `maximal_polygons :: Vector{Vector{RationalPolygon{T}}}`,
- `all_polygons :: Vector{Vector{RationalPolygon{T}}}`,
- `total_count :: Int`.

"""
mutable struct InMemoryCastryckStorage{T <: Integer} <: CastryckStorage{T}
    maximum_genus :: Int
    maximal_polygons :: Vector{Vector{RationalPolygon{T}}}
    all_polygons :: Vector{Vector{RationalPolygon{T}}}
    total_count :: Int

    function InMemoryCastryckStorage{T}(maximum_genus :: Int = 100) where {T <: Integer}
        maximal_polygons = Vector{RationalPolygon{T}}[]
        push!(maximal_polygons, classify_maximal_polygons_genus_one(1))
        for i = 2 : maximum_genus
            push!(maximal_polygons, RationalPolygon{T}[])
        end
        all_polygons = Vector{RationalPolygon{T}}[]
        return new{T}(maximum_genus, maximal_polygons, all_polygons, 0)
    end

end

last_completed_genus(st :: InMemoryCastryckStorage{T}) where {T <: Integer} =
length(st.all_polygons)


@doc raw"""
    classify_next_genus(st :: InMemoryCastryckStorage{T}; logging :: Bool = false) where {T <: Integer}

Perform a single step in Castryck's classification of lattice polygons by
number of interior lattice points, using in-memory storage. Returns a tuple
where the first entry is the number of lattice polygons obtained and the second
number is the number of maximal lattice polygons.

# Example

Perform two steps in Castryck's classification. The result tells us that there
are 45 lattice polygons with exactly two interior lattice points, four of which
are maximal.


```jldoctest
julia> st = InMemoryCastryckStorage{Int}();

julia> classify_next_genus(st)
(16, 3)

julia> classify_next_genus(st)
(45, 4)
```

"""
function classify_next_genus(st :: InMemoryCastryckStorage{T}; logging :: Bool = false) where {T <: Integer}

    i = last_completed_genus(st) + 1

    Ps_max = st.maximal_polygons[i]

    append!(Ps_max, classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(i,T))
    append!(Ps_max, classify_maximal_lattice_polygons_with_collinear_interior_points(i,T))

    logging && @info "[i = $i]. Got $(length(Ps_max)) maximal polygons."

    new_polygons = subpolygons(Ps_max;
        only_equal_number_of_interior_lattice_points = true,
        use_affine_normal_form = true)

    push!(st.all_polygons, new_polygons)
    st.total_count += length(new_polygons)

    logging && @info "[i = $i]. Subpolygons complete. Num of polygons: $(length(new_polygons))"

    moved_out_polygons = Dict{Int, Vector{RationalPolygon{T}}}[]
    for i = 1 : Threads.nthreads()
        push!(moved_out_polygons, Dict{Int, Vector{RationalPolygon{T}}}())
    end

    Threads.@threads for P ∈ new_polygons
        _is_internal_quick(P) || continue
        n = number_of_lattice_points(P)
        n ≤ st.maximum_genus || continue
        Q = move_out_edges(P)
        rationality(Q) == 1 || continue

        tid = Threads.threadid()
        if !haskey(moved_out_polygons[tid], n)
            moved_out_polygons[tid][n] = RationalPolygon{T}[]
        end
        push!(moved_out_polygons[tid][n], Q)
    end

    total_moved_out_polygons = mergewith(append!, moved_out_polygons...)

    for (n,Ps) ∈ total_moved_out_polygons
        append!(st.maximal_polygons[n], Ps)
    end

    logging && @info "[i = $i]. Moving out complete. New maximal polygons: $(sum(length.(values(total_moved_out_polygons))))"

    return length(new_polygons), length(Ps_max)

end


@doc raw"""
    struct HDFCastryckStoragePreferences{T <: Integer}

A struct holding preferences for Castryck's classification using the HDF5 file format. It has four fields:

- `swmr :: Bool`: Whether to use single-reader-multiple-writer mode for HDF5.
    Defaults to `true`.
- `maximum_genus :: Int`: An upper bound for the maximal number of interior
    lattice points to which the classification should be run. Defaults to `100`.
- `maximum_number_of_vertices :: Int`: An upper bound for the maximal number of
    vertices to be expected in the classification. This has to be set since every
    HDF5 file generated will have a dataset "numbers\_of\_polygons" storing the
    number of polygons for each number of vertices and the size of this dataset
    needs to be set beforehand. Defaults to `100`, which should be more than enough for any
    feasable computation.
- `block_size :: Int`: How many polygons should be read into memory at once
    during the computation of subpolygons and the moving-out process. Defaults to
    `10^6`.


"""
struct HDFCastryckStoragePreferences{T <: Integer}
    swmr :: Bool
    maximum_genus :: Int
    maximum_number_of_vertices :: Int
    block_size :: Int
    
    HDFCastryckStoragePreferences{T}(
        swmr :: Bool = true,
        maximum_genus :: Int = 100,
        maximum_number_of_vertices :: Int = 100,
        block_size :: Int = 10^6) where {T <: Integer} =
    new{T}(swmr, maximum_genus, maximum_number_of_vertices, block_size)

end


@doc raw"""
    mutable struct HDFCastryckStorage{T <: Integer} <: CastryckStorage{T}

A struct for managing classification results of Castryck's classification of lattice polygons using the HDF5 file format. It has the following fields:

- `preferences :: HDFCastryckStoragePreferences{T}`
- `directory_path :: String`: The directory where the HDF5 files will be generated.
- `last_completed_genus :: Int`: The last completed step of the classification. Initially, this will be `0`.
- `total_count :: Int`

"""
mutable struct HDFCastryckStorage{T <: Integer} <: CastryckStorage{T}
    preferences :: HDFCastryckStoragePreferences{T}
    directory_path :: String
    last_completed_genus :: Int
    total_count :: Int

    function HDFCastryckStorage{T}(
        preferences :: HDFCastryckStoragePreferences{T},
        directory_path :: String) where {T <: Integer}

        isdir(directory_path) || error("$directory_path is not a directory")
        mkdir(joinpath(directory_path, "maximal"))
        mkdir(joinpath(directory_path, "all"))

        f = h5open(joinpath(directory_path, "maximal/i1.h5"), "cw"; swmr = preferences.swmr)

        Ps3 = [RationalPolygon(SMatrix{2,3,T}(1,0,1,2,-3,-4), 1), RationalPolygon(SMatrix{2,3,T}(1,0,1,3,-2,-3), 1)]
        Ps4 = [RationalPolygon(SMatrix{2,4,T}(1,0,1,2,-1,0,-1,-2), 1)]
        write_polygon_dataset(f, "n3", Ps3)
        write_polygon_dataset(f, "n4", Ps4)

        f["numbers_of_polygons"] = zeros(Int, preferences.maximum_number_of_vertices) 
        f["numbers_of_polygons"][3] = 2
        f["numbers_of_polygons"][4] = 1

        close(f)

        for i = 2 : preferences.maximum_number_of_vertices
            f = h5open(joinpath(directory_path, "maximal/i$i.h5"), "cw"; swmr = preferences.swmr)
            f["numbers_of_polygons"] = zeros(Int, preferences.maximum_number_of_vertices) 
            close(f)
        end

        new{T}(preferences, directory_path, 0, 0)

    end

    HDFCastryckStorage{T}(
        directory_path :: String;
        swmr :: Bool = true,
        maximum_genus :: Int = 100,
        maximum_number_of_vertices :: Int = 100,
        block_size :: Int = 10^6) where {T <: Integer} =
    HDFCastryckStorage{T}(HDFCastryckStoragePreferences{T}(swmr, maximum_genus, maximum_number_of_vertices, block_size), directory_path)

end

last_completed_genus(st :: HDFCastryckStorage{T}) where {T <: Integer} =
st.last_completed_genus


@doc raw"""
    classify_next_genus(st :: HDFCastryckStorage{T}; logging :: Bool = false) where {T <: Integer}

Perform a single step in Castryck's classification of lattice polygons using
on-disk storage with HDF5.

"""
function classify_next_genus(st :: HDFCastryckStorage{T}; logging :: Bool = false) where {T <: Integer}

    i = last_completed_genus(st) + 1

    Ps_max = RationalPolygon{T}[]
    f = h5open(joinpath(st.directory_path, "maximal/i$i.h5"), "r"; swmr = st.preferences.swmr)
    for n ∈ keys(f)
        n != "numbers_of_polygons" || continue
        append!(Ps_max, read_polygon_dataset(one(T), f, n))
    end

    append!(Ps_max, classify_maximal_lattice_polygons_with_two_dimensional_empty_fine_interior(i,T))
    append!(Ps_max, classify_maximal_lattice_polygons_with_collinear_interior_points(i,T))

    logging && @info "[i = $i]. Got $(length(Ps_max)) maximal polygons."

    subpolygons_storage = HDFSubpolygonStorage{T}(
        Ps_max,
        joinpath(st.directory_path, "all/i$i.h5");
        use_affine_normal_form = true,
        only_equal_number_of_interior_lattice_points = true,
        swmr = st.preferences.swmr,
        block_size = st.preferences.block_size)

    subpolygons(subpolygons_storage)
    st.total_count += subpolygons_storage.total_count

    @info "[i = $i]. Subpolygons complete. Num of polygons: $(subpolygons_storage.total_count)"

    moved_out_polygons = Dict{Tuple{Int,Int}, Vector{RationalPolygon{T}}}[]
    for i = 1 : Threads.nthreads()
        push!(moved_out_polygons, Dict{Tuple{Int,Int}, Vector{RationalPolygon{T}}}())
    end

    f = h5open(joinpath(st.directory_path, "all/i$i.h5"), "r"; swmr = st.preferences.swmr)
    block_size = st.preferences.block_size
    for a_string ∈ keys(f)
        startswith(a_string, "a") || continue
        a = parse(Int, a_string[2:end])
        for n_string ∈ keys(f[a_string])
            N = length(dataspace(f[a_string][n_string]))
            number_of_blocks = N ÷ block_size + 1

            for b = 1 : number_of_blocks
                I = (b-1)*block_size + 1 : min(b*block_size, N)
                new_count = 0
                n = parse(Int, n_string[2:end])
                Ps = read_polygon_dataset(one(T), f[a_string], n_string, I) 

                Threads.@threads for P ∈ Ps
                    _is_internal_quick(P) || continue
                    iQ = number_of_lattice_points(P)
                    iQ ≤ st.preferences.maximum_genus || continue
                    Q = move_out_edges(P)
                    rationality(Q) == 1 || continue
                    nQ = number_of_vertices(Q)

                    tid = Threads.threadid()
                    if !haskey(moved_out_polygons[tid], (iQ,nQ))
                        moved_out_polygons[tid][(iQ,nQ)] = RationalPolygon{T}[]
                    end
                    push!(moved_out_polygons[tid][(iQ,nQ)], Q)
                end
            end
        end
    end

    total_moved_out_polygons = mergewith(append!, moved_out_polygons...)
    new_maximal_polygons_count = sum(length.(values(total_moved_out_polygons)))

    logging && @info "[i = $i]. Moving out complete. New maximal polygons: $new_maximal_polygons_count"

    for iQ ∈ unique(first.(keys(total_moved_out_polygons)))
        f = h5open(joinpath(st.directory_path, "maximal/i$iQ.h5"), "r+"; swmr = st.preferences.swmr)
        for nQ ∈ unique([k[2] for k ∈ keys(total_moved_out_polygons) if k[1] == iQ])
            Qs = convert(Vector{RationalPolygon{T,nQ,2nQ}}, total_moved_out_polygons[(iQ,nQ)])
            write_polygon_dataset(f, "n$nQ", Qs)
            f["numbers_of_polygons"][nQ] += length(Qs)
        end
        close(f)
    end

    st.last_completed_genus = i

    return subpolygons_storage.total_count, length(Ps_max)

end


@doc raw"""
    classify_lattice_polygons_by_genus(st :: CastryckStorage{T}, max_genus :: Int; logging :: Bool = false) where {T <: Integer}

Run Castryck's classification of lattice polygons by number of interior lattice
points, up to `max_genus`. The classification is multithreaded, so make sure
julia has access to a good number of threads for maximum performance (i.e.
`Threads.nthreads()` is greater than one).

# Example

Reproduce Castryck's classification in memory, see Table 1 of [Cas12](@cite) or
A322343 on OEIS. This should not take longer than a few minutes on modern
hardware.

```julia
julia> st = InMemoryCastryckStorage{Int}();

julia> classify_lattice_polygons_by_genus(st, 30; logging=true)
[ Info: [i = 1]. Got 3 maximal polygons.
[ Info: [i = 1]. Subpolygons complete. Num of polygons: 16
[ Info: [i = 1]. Moving out complete. New maximal polygons: 16
[ Info: [i = 2]. Got 4 maximal polygons.
[ Info: [i = 2]. Subpolygons complete. Num of polygons: 45
[ Info: [i = 2]. Moving out complete. New maximal polygons: 22
[ Info: [i = 3]. Got 6 maximal polygons.
[ Info: [i = 3]. Subpolygons complete. Num of polygons: 120
[ Info: [i = 3]. Moving out complete. New maximal polygons: 63
[ Info: [i = 4]. Got 9 maximal polygons.
[ Info: [i = 4]. Subpolygons complete. Num of polygons: 211
[ Info: [i = 4]. Moving out complete. New maximal polygons: 78
[ Info: [i = 5]. Got 11 maximal polygons.
[ Info: [i = 5]. Subpolygons complete. Num of polygons: 403
[ Info: [i = 5]. Moving out complete. New maximal polygons: 122
[ Info: [i = 6]. Got 13 maximal polygons.
[ Info: [i = 6]. Subpolygons complete. Num of polygons: 714
[ Info: [i = 6]. Moving out complete. New maximal polygons: 192
[ Info: [i = 7]. Got 16 maximal polygons.
[ Info: [i = 7]. Subpolygons complete. Num of polygons: 1023
[ Info: [i = 7]. Moving out complete. New maximal polygons: 239
[ Info: [i = 8]. Got 21 maximal polygons.
[ Info: [i = 8]. Subpolygons complete. Num of polygons: 1830
[ Info: [i = 8]. Moving out complete. New maximal polygons: 316
[ Info: [i = 9]. Got 27 maximal polygons.
[ Info: [i = 9]. Subpolygons complete. Num of polygons: 2700
[ Info: [i = 9]. Moving out complete. New maximal polygons: 508
[ Info: [i = 10]. Got 33 maximal polygons.
[ Info: [i = 10]. Subpolygons complete. Num of polygons: 3659
[ Info: [i = 10]. Moving out complete. New maximal polygons: 509
[ Info: [i = 11]. Got 38 maximal polygons.
[ Info: [i = 11]. Subpolygons complete. Num of polygons: 6125
[ Info: [i = 11]. Moving out complete. New maximal polygons: 700
[ Info: [i = 12]. Got 51 maximal polygons.
[ Info: [i = 12]. Subpolygons complete. Num of polygons: 8101
[ Info: [i = 12]. Moving out complete. New maximal polygons: 1044
[ Info: [i = 13]. Got 61 maximal polygons.
[ Info: [i = 13]. Subpolygons complete. Num of polygons: 11027
[ Info: [i = 13]. Moving out complete. New maximal polygons: 1113
[ Info: [i = 14]. Got 76 maximal polygons.
[ Info: [i = 14]. Subpolygons complete. Num of polygons: 17280
[ Info: [i = 14]. Moving out complete. New maximal polygons: 1429
[ Info: [i = 15]. Got 86 maximal polygons.
[ Info: [i = 15]. Subpolygons complete. Num of polygons: 21499
[ Info: [i = 15]. Moving out complete. New maximal polygons: 2052
[ Info: [i = 16]. Got 113 maximal polygons.
[ Info: [i = 16]. Subpolygons complete. Num of polygons: 28689
[ Info: [i = 16]. Moving out complete. New maximal polygons: 1962
[ Info: [i = 17]. Got 129 maximal polygons.
[ Info: [i = 17]. Subpolygons complete. Num of polygons: 43012
[ Info: [i = 17]. Moving out complete. New maximal polygons: 2651
[ Info: [i = 18]. Got 166 maximal polygons.
[ Info: [i = 18]. Subpolygons complete. Num of polygons: 52736
[ Info: [i = 18]. Moving out complete. New maximal polygons: 3543
[ Info: [i = 19]. Got 200 maximal polygons.
[ Info: [i = 19]. Subpolygons complete. Num of polygons: 68557
[ Info: [i = 19]. Moving out complete. New maximal polygons: 3638
[ Info: [i = 20]. Got 240 maximal polygons.
[ Info: [i = 20]. Subpolygons complete. Num of polygons: 97733
[ Info: [i = 20]. Moving out complete. New maximal polygons: 4594
[ Info: [i = 21]. Got 281 maximal polygons.
[ Info: [i = 21]. Subpolygons complete. Num of polygons: 117776
[ Info: [i = 21]. Moving out complete. New maximal polygons: 5996
[ Info: [i = 22]. Got 352 maximal polygons.
[ Info: [i = 22]. Subpolygons complete. Num of polygons: 152344
[ Info: [i = 22]. Moving out complete. New maximal polygons: 6364
[ Info: [i = 23]. Got 403 maximal polygons.
[ Info: [i = 23]. Subpolygons complete. Num of polygons: 209409
[ Info: [i = 23]. Moving out complete. New maximal polygons: 7922
[ Info: [i = 24]. Got 506 maximal polygons.
[ Info: [i = 24]. Subpolygons complete. Num of polygons: 248983
[ Info: [i = 24]. Moving out complete. New maximal polygons: 9692
[ Info: [i = 25]. Got 584 maximal polygons.
[ Info: [i = 25]. Subpolygons complete. Num of polygons: 319957
[ Info: [i = 25]. Moving out complete. New maximal polygons: 10208
[ Info: [i = 26]. Got 708 maximal polygons.
[ Info: [i = 26]. Subpolygons complete. Num of polygons: 420714
[ Info: [i = 26]. Moving out complete. New maximal polygons: 12727
[ Info: [i = 27]. Got 821 maximal polygons.
[ Info: [i = 27]. Subpolygons complete. Num of polygons: 497676
[ Info: [i = 27]. Moving out complete. New maximal polygons: 15431
[ Info: [i = 28]. Got 995 maximal polygons.
[ Info: [i = 28]. Subpolygons complete. Num of polygons: 641229
[ Info: [i = 28]. Moving out complete. New maximal polygons: 15918
[ Info: [i = 29]. Got 1121 maximal polygons.
[ Info: [i = 29]. Subpolygons complete. Num of polygons: 813814
[ Info: [i = 29]. Moving out complete. New maximal polygons: 20354
[ Info: [i = 30]. Got 1352 maximal polygons.
[ Info: [i = 30]. Subpolygons complete. Num of polygons: 957001
[ Info: [i = 30]. Moving out complete. New maximal polygons: 23873
```

# Example

Reproduce Castryck's classification and storing the output to HDF5 files

```julia
julia> st = HDFCastryckStorage{Int}("/tmp");

julia> classify_lattice_polygons_by_genus(st, 30);
```

This will create two directories "all" and "maximal" in the target directory
and populate them with HDF5 files "i1.h5, i2.h5, ..." containing the polygons.
The files in "maximal" contain datasets ordered by number of vertices. The
files in "all" contain groups ordered by area, which contain datasets ordered
by number of vertices. Every h5 file additionally contains a dataset
"numbers\_of\_polygons".

```julia
julia> using HDF5

julia> f = h5open("/tmp/all/i30.h5", "r");

julia> Ps = read_polygon_dataset(1, f, "a69/n8"); # Read in all octagons of genus 30 with normalized volume 69

julia> all(P -> number_of_interior_lattice_points(P) == 30, Ps)
true

julia> A = read_dataset(f, "numbers_of_polygons"); # The total number of lattice polygons with 30 interior lattice points

julia> sum(A) # The total number of lattice polygons with 30 interior lattice points.
957001
```

"""
function classify_lattice_polygons_by_genus(st :: CastryckStorage{T}, max_genus :: Int; logging :: Bool = false) where {T <: Integer}
    for i = last_completed_genus(st) + 1 : max_genus
        new_total_count, new_max_count = classify_next_genus(st; logging)
    end
end


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
