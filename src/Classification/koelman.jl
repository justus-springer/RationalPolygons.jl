@doc raw"""
    height_one_points(P :: RationalPolygon)

Given a lattice polygon `P`, return all lattice points that have lattice height
one with respect to some edge of `P`. Equivalently, return the set of boundary
lattice points of `move_out_edges(P)`.

"""
function height_one_points(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    V = vertex_matrix(P)
    Hs = affine_halfplanes(P)
    lines = [line(Hs[i] - 1) for i = 1 : N]
    res = LatticePoint{T}[]
    for i = 1 : N
        p = intersection_point(lines[mod(i-1,1:N)], lines[i])
        q = intersection_point(lines[mod(i+1,1:N)], lines[i])
        nv = normal_vector(lines[i])
        _, a, b = gcdx(-nv[1],-nv[2])
        x0 = LatticePoint{T}(V[1,i] + a, V[2,i] + b)
        for b ∈ integral_points_on_line_segment_with_given_integral_point(p,q,x0)
            all(H -> distance(b,H) ≥ -1, Hs) || continue
            b ∉ res || continue
            push!(res, b)
        end
    end
    return res
end


@doc raw"""
    single_point_extensions(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

Return all lattice polygons that can be obtained by adding a single height one
point to a polygon of `Ps`, up to affine equivalence.

"""
function single_point_extensions(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    out_dicts = Vector{Dict{Int, Set{<:RationalPolygon{T}}}}(undef, Threads.nthreads())
    for i = 1 : Threads.nthreads()
        out_dicts[i] = Dict{Int, Set{<:RationalPolygon{T}}}()
    end

    Threads.@threads for P ∈ Ps
        tid = Threads.threadid()
        N = number_of_vertices(P)
        bs = height_one_points(P)
        V = vertex_matrix(P)
        Hs = affine_halfplanes(P)
        for b ∈ bs
            contained_in_last_halfplane = contains_in_interior(b, Hs[N])
            local entry_point :: Int
            local exit_point :: Int
            for i = 1 : N
                contained_in_this_halfplane = contains_in_interior(b, Hs[i])
                if contained_in_last_halfplane && !contained_in_this_halfplane
                    exit_point = i
                elseif contained_in_this_halfplane && !contained_in_last_halfplane
                    entry_point = i
                end
                contained_in_last_halfplane = contained_in_this_halfplane
            end
            exit_point < entry_point && (exit_point += N)
            new_vertices = LatticePoint{T}[b]
            for i = entry_point : exit_point
                push!(new_vertices, V[:,mod(i,1:N)])
            end
            Q = RationalPolygon(new_vertices, one(T))
            (Q,is_special) = lattice_affine_normal_form_with_is_first_index_special(Q)
            is_special || continue
            
            M = number_of_vertices(Q)
            !haskey(out_dicts[tid], M) && (out_dicts[tid][M] = Set{RationalPolygon{T,M,2M}}())
            push!(out_dicts[tid][M], Q)
        end
    end

    for i = 2 : length(out_dicts)
        mergewith!(union!, out_dicts[1], out_dicts[i])
        out_dicts[i] = Dict{Int,Set{<:RationalPolygon{T}}}()
    end

    return out_dicts[1]
end


@doc raw"""
    abstract type KoelmanStorage{T <: Integer} end

Abstract supertype of `InMemoryKoelmanStorage` and `HDFKoelmanStorage`. Both
implement `classify_next_number_of_lattice_points`, which performs a single
step in Koelman's classification of lattice polygons.

"""
abstract type KoelmanStorage{T <: Integer} end


@doc raw"""
    mutable struct InMemoryKoelmanStorage{T <: Integer} <: KoelmanStorage{T}

A struct holding classification results of Koelman's classification of lattice
polygons by number of lattice points. It has two fields `polygons ::
Vector{Vector{RationalPolygon{T}}` and `total_count :: Int`.


"""
mutable struct InMemoryKoelmanStorage{T <: Integer} <: KoelmanStorage{T}
    polygons :: Vector{Vector{RationalPolygon{T}}}
    total_count :: Int

    function InMemoryKoelmanStorage{T}() where {T <: Integer}
        polygons = [RationalPolygon{T}[], RationalPolygon{T}[], [RationalPolygon(SMatrix{2,3,T,6}(0,0,1,0,0,1), one(T))]]
        total_count = 1
        return new{T}(polygons, total_count)
    end
end

last_completed_number_of_lattice_points(st :: InMemoryKoelmanStorage{T}) where {T <: Integer} = length(st.polygons)


@doc raw"""
    classify_next_number_of_lattice_points(st :: InMemoryKoelmanStorage{T}) where {T <: Integer}

Perform a single step in Koelman's classification of lattice polygons using
in-memory storage.

# Example

Perform a single step of Koelmans classification using `Int64`. The result
tells us that there are three lattice polygons with exactly four lattice
points.

```jldoctest
julia> st = InMemoryKoelmanStorage{Int64}();

julia> classify_next_number_of_lattice_points(st)
3

julia> st.polygons[4]
3-element Vector{RationalPolygon{Int64}}:
 Rational polygon of rationality 1 with 3 vertices.
 Rational polygon of rationality 1 with 4 vertices.
 Rational polygon of rationality 1 with 3 vertices.

```

"""
function classify_next_number_of_lattice_points(st :: InMemoryKoelmanStorage{T}; logging :: Bool = false) where {T <: Integer}
    Ps = last(st.polygons)
    new_Ps = collect(union(values(single_point_extensions(Ps))...))
    push!(st.polygons, new_Ps)
    new_count = length(new_Ps)
    st.total_count += new_count
    return new_count
end



classify_polygons_by_number_of_lattice_points(max_number_of_lattice_points :: Int, T :: Type{<:Integer} = Int; logging :: Bool = false) =
classify_polygons_by_number_of_lattice_points(InMemoryKoelmanStorage{T}(), max_number_of_lattice_points; logging)


@doc raw"""
    struct HDFKoelmanStoragePreferences{T <: Integer}   

A struct holding preferences for Koelman's classification using the HDF5 file format. It has four fields:

- `swmr :: Bool`: Whether to use single-reader-multiple-writer mode for HDF5.
    Defaults to `true`.
- `maximum_number_of_vertices :: Int`: An upper bound for the maximal number of
    vertices to be expected in the classification. This has to be set since every
    HDF5 file generated will have a dataset "numbers\_of\_polygons" storing the
    number of polygons for each number of vertices and the size of this dataset
    needs to be set beforehand. Defaults to `100`, which should be more than enough for any
    feasable computation.
- `block_size :: Int`: How many polygons should be read into memory at once
    during the extension process. Defaults to `10^6`.

"""
struct HDFKoelmanStoragePreferences{T <: Integer}
    swmr :: Bool
    maximum_number_of_vertices :: Int
    block_size :: Int
    
    HDFKoelmanStoragePreferences{T}(
        swmr :: Bool = true,
        maximum_number_of_vertices :: Int = 100,
        block_size :: Int = 10^6) where {T <: Integer} =
    new{T}(swmr, maximum_number_of_vertices, block_size)

end


@doc raw"""
    mutable struct HDFKoelmanStorage{T <: Integer} <: KoelmanStorage{T}

A struct for managing classification results of Koelman's classification of lattice polygons using the HDF5 file format. It has the following fields:

- `preferences :: HDFKoelmanStoragePreferences{T}`
- `directory_path :: String`: The directory where the HDF5 files will be generated.
- `last_completed_number_of_lattice_points :: Int`: The last completed step of the classification. Initially, this will be `3`.
- `total_count :: Int`

"""
mutable struct HDFKoelmanStorage{T <: Integer} <: KoelmanStorage{T}
    preferences :: HDFKoelmanStoragePreferences{T}
    directory_path :: String
    last_completed_number_of_lattice_points :: Int
    total_count :: Int

    function HDFKoelmanStorage{T}(
        preferences :: HDFKoelmanStoragePreferences{T},
        directory_path :: String) where {T <: Integer}

        f = h5open(joinpath(directory_path, "l3.h5"), "cw"; swmr = preferences.swmr)

        P = RationalPolygon(SMatrix{2,3,T,6}(0,0,1,0,0,1), 1)
        write_polygon_dataset(f, "n3", [P])

        f["numbers_of_polygons"] = zeros(Int, preferences.maximum_number_of_vertices) 
        f["numbers_of_polygons"][3] = 1

        close(f)

        new{T}(preferences, directory_path, 3, 1)

    end

    HDFKoelmanStorage{T}(directory_path :: String;
                         swmr :: Bool = true,
                         maximum_number_of_vertices :: Int = 100,
                         block_size :: Int = 10^6) where {T <: Integer} =
    HDFKoelmanStorage{T}(HDFKoelmanStoragePreferences{T}(swmr, maximum_number_of_vertices, block_size), directory_path)

end

last_completed_number_of_lattice_points(st :: HDFKoelmanStorage{T}) where {T <: Integer} =
st.last_completed_number_of_lattice_points


@doc raw"""
    classify_next_number_of_lattice_points(st :: HDFKoelmanStorage{T}; logging :: Bool = false) where {T <: Integer}

Perform a single step in Koelman's classification of lattice polygons using
on-disk storage with HDF5.

"""
function classify_next_number_of_lattice_points(st :: HDFKoelmanStorage{T}; logging :: Bool = false) where {T <: Integer}

    l = st.last_completed_number_of_lattice_points
    block_size = st.preferences.block_size

    last_file = h5open(joinpath(st.directory_path, "l$l.h5"), "r"; swmr = st.preferences.swmr)
    current_file = h5open(joinpath(st.directory_path, "l$(l+1).h5"), "cw"; swmr = st.preferences.swmr)
    current_file["numbers_of_polygons"] = zeros(Int, st.preferences.maximum_number_of_vertices) 

    total_new_count = 0

    for n_string ∈ keys(last_file)
        n_string != "numbers_of_polygons" || continue

        N = length(dataspace(last_file[n_string]))
        number_of_blocks = N ÷ block_size + 1

        for b = 1 : number_of_blocks
            elapsed_time = @elapsed begin

                I = (b-1)*block_size + 1 : min(b*block_size, N)
                new_count = 0
                n = parse(Int, n_string[2:end])
                Ps = read_polygon_dataset(one(T), last_file, n_string, I) 

                logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Polygons to extend: $(length(Ps))"

                new_Ps = single_point_extensions(Ps)

                logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Extension complete. New polygons: $(sum(length.(values(new_Ps))))"

                for (m, Qs) ∈ new_Ps
                    write_polygon_dataset(current_file, "n$m", collect(Qs))
                    new_count += length(Qs)
                    total_new_count += length(Qs)
                    current_file["numbers_of_polygons"][m] += length(Qs)
                    st.total_count += length(Qs)
                end
                logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Writeout complete. Running total: $(st.total_count)"
            end
            pps = round(new_count / elapsed_time, digits=2)
            elapsed_time = round(elapsed_time, digits=2)
            logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Last step took $(elapsed_time) seconds. Averaging $pps polygons/s"
        end
    end

    st.last_completed_number_of_lattice_points = l+1

    close(last_file)
    close(current_file)

    return total_new_count

end


@doc raw"""
    classify_polygons_by_number_of_lattice_points(st :: KoelmanStorage{T}, max_number_of_lattice_points :: Int; logging :: Bool = false) where {T <: Integer}

Run Koelman's classification of lattice polygons by number of lattice points,
up to `max_number_of_lattice_points`. The classification is multithreaded, so
make sure julia has access to a good number of threads for maximum performance
(i.e. `Threads.nthreads()` is greater than one).

# Example

Reproduce Koelman's original classification in memory, see Table 4.4.3 of
[Koe91](@cite) or A371917 on OEIS. This should not take longer than a few
minutes on modern hardware.

```julia
julia> st = InMemoryKoelmanStorage{Int}();

julia> classify_polygons_by_number_of_lattice_points(st, 42; logging=true);
[ Info: [l = 4]. New polygons: 3. Total: 4
[ Info: [l = 5]. New polygons: 6. Total: 10
[ Info: [l = 6]. New polygons: 13. Total: 23
[ Info: [l = 7]. New polygons: 21. Total: 44
[ Info: [l = 8]. New polygons: 41. Total: 85
[ Info: [l = 9]. New polygons: 67. Total: 152
[ Info: [l = 10]. New polygons: 111. Total: 263
[ Info: [l = 11]. New polygons: 175. Total: 438
[ Info: [l = 12]. New polygons: 286. Total: 724
[ Info: [l = 13]. New polygons: 419. Total: 1143
[ Info: [l = 14]. New polygons: 643. Total: 1786
[ Info: [l = 15]. New polygons: 938. Total: 2724
[ Info: [l = 16]. New polygons: 1370. Total: 4094
[ Info: [l = 17]. New polygons: 1939. Total: 6033
[ Info: [l = 18]. New polygons: 2779. Total: 8812
[ Info: [l = 19]. New polygons: 3819. Total: 12631
[ Info: [l = 20]. New polygons: 5293. Total: 17924
[ Info: [l = 21]. New polygons: 7191. Total: 25115
[ Info: [l = 22]. New polygons: 9752. Total: 34867
[ Info: [l = 23]. New polygons: 12991. Total: 47858
[ Info: [l = 24]. New polygons: 17321. Total: 65179
[ Info: [l = 25]. New polygons: 22641. Total: 87820
[ Info: [l = 26]. New polygons: 29687. Total: 117507
[ Info: [l = 27]. New polygons: 38533. Total: 156040
[ Info: [l = 28]. New polygons: 49796. Total: 205836
[ Info: [l = 29]. New polygons: 63621. Total: 269457
[ Info: [l = 30]. New polygons: 81300. Total: 350757
[ Info: [l = 31]. New polygons: 102807. Total: 453564
[ Info: [l = 32]. New polygons: 129787. Total: 583351
[ Info: [l = 33]. New polygons: 162833. Total: 746184
[ Info: [l = 34]. New polygons: 203642. Total: 949826
[ Info: [l = 35]. New polygons: 252898. Total: 1202724
[ Info: [l = 36]. New polygons: 313666. Total: 1516390
[ Info: [l = 37]. New polygons: 386601. Total: 1902991
[ Info: [l = 38]. New polygons: 475540. Total: 2378531
[ Info: [l = 39]. New polygons: 582216. Total: 2960747
[ Info: [l = 40]. New polygons: 710688. Total: 3671435
[ Info: [l = 41]. New polygons: 863552. Total: 4534987
[ Info: [l = 42]. New polygons: 1048176. Total: 5583163
```

# Example

Reproduce Koelman's classification and storing the output to HDF5 files.

```julia
julia> st = HDFKoelmanStorage{Int}("/tmp");

julia> classify_polygons_by_number_of_lattice_points(st, 42);
```

The result will be HDF5 files named "l1.h5, l2.h5, ...", each containing one
dataset of polygons for fixed number of vertices as well as a dataset
"numbers\_of\_polygons" holding their numbers.

```julia
julia> using HDF5

julia> f = h5open("/tmp/test/l42.h5", "r")
🗂️ HDF5.File: (read-only) /tmp/test/l42.h5
├─ 🔢 n10
├─ 🔢 n11
├─ 🔢 n12
├─ 🔢 n3
├─ 🔢 n4
├─ 🔢 n5
├─ 🔢 n6
├─ 🔢 n7
├─ 🔢 n8
├─ 🔢 n9
└─ 🔢 numbers_of_polygons

julia> A = read_dataset(f, "numbers_of_polygons");

julia> sum(A) # the number of lattice polygons with 42 lattice points
1048176

julia> Ps = read_polygon_dataset(1, f, "n5"); # Read in all pentagons with 42 lattice points.

julia> all(P -> number_of_lattice_points(P) == 42, Ps)
true
```

"""
function classify_polygons_by_number_of_lattice_points(st :: KoelmanStorage{T}, max_number_of_lattice_points :: Int; logging :: Bool = false) where {T <: Integer}
    for l = last_completed_number_of_lattice_points(st) + 1 : max_number_of_lattice_points
        new_count = classify_next_number_of_lattice_points(st; logging)
        logging && @info "[l = $l]. New polygons: $new_count. Total: $(st.total_count)"
    end
end
