
@doc raw"""
    classify_lattice_triangles_by_picard_index(p :: T) where {T <: Integer}

Return all lattice triangles with Picard index `p`.

# Example

There are two lattice triangles with Picard index 6:

```jldoctest
julia> classify_lattice_triangles_by_picard_index(6)
Set{RationalPolygon{Int64, 3, 6}} with 2 elements:
  Rational polygon of rationality 1 with 3 vertices.
  Rational polygon of rationality 1 with 3 vertices.
```

"""
function classify_lattice_triangles_by_picard_index(p :: T) where {T <: Integer}
    res = Set{RationalPolygon{T,3,6}}()

    for n = 1 : isqrt(p)
        q0,r0 = divrem(p, n^2)
        r0 == 0 || continue
        for w0 = 1 : isqrt(q0)
            q1,r1 = divrem(q0, w0)
            r1 == 0 || continue
            gcd(w0, q0 ÷ w0) == 1 || continue
            for w1 = w0 : isqrt(q1)
                w2,r2 = divrem(q1, w1)
                r2 == 0 || continue
                gcd(w1,w2) == 1 || continue
                for x = 0 : n*w0-1
                    (y,r) = divrem(w2 + x*w1, w0)
                    r == 0 || continue
                    gcd(x, n*w0) == 1 || continue
                    gcd(y, n*w1) == 1 || continue
                    P = RationalPolygon(SMatrix{2,3,T,6}(1,0,x,n*w0,-y,-n*w1), 1)
                    push!(res, unimodular_normal_form(P))
                end
            end
        end
    end
    return res
end


@doc raw"""
    abstract type PicardIndexStorage{T <: Integer} end

Abstract supertype of `InMemoryPicardIndexStorage` and `HDFPicardIndexStorage`.

"""
abstract type PicardIndexStorage{T <: Integer} end


@doc raw"""
    mutable struct InMemoryPicardIndexStorage{T <: Integer} <: PicardIndexStorage{T}

A struct holding classification results of Springer's classification of lattice
triangles by picard index.

"""
mutable struct InMemoryPicardIndexStorage{T <: Integer} <: PicardIndexStorage{T}
    polygons :: Vector{Vector{RationalPolygon{T,3,6}}}
    last_completed_picard_index :: T

    InMemoryPicardIndexStorage{T}() where {T <: Integer} = new{T}(Vector{RationalPolygon{T,3,6}}[], 0)

end


@doc raw"""
    classify_lattice_triangles_by_picard_index(st :: InMemoryPicardIndexStorage{T}, max_picard_index :: T) where {T <: Integer}

Perform Springer's classification of lattice triangles up go
`max_picard_index`, storing the results in memory.

# Example

Reproduce Springer's original classification up to picard index 10000, see
Theorem 1.2 of [Spr24](@cite).

```jldoctest
julia> st = InMemoryPicardIndexStorage{Int}()
InMemoryPicardIndexStorage{Int64}(Vector{RationalPolygon{Int64, 3, 6}}[], 0)

julia> classify_lattice_triangles_by_picard_index(st, 10000);

julia> sum(length.(st.polygons))
68053
```

"""
function classify_lattice_triangles_by_picard_index(st :: InMemoryPicardIndexStorage{T}, max_picard_index :: T) where {T <: Integer}
    dicts = Dict{T, Set{RationalPolygon{T,3,6}}}[]
    for t = 1 : Threads.nthreads()
        push!(dicts, Dict{T, Set{RationalPolygon{T,3,6}}}())
    end

    p0 = st.last_completed_picard_index + 1
    Threads.@threads for p = p0 : max_picard_index
        t = Threads.threadid()
        dicts[t][p] = classify_lattice_triangles_by_picard_index(p)
    end

    merged_dict = merge!(dicts...)
    for p = p0 : max_picard_index
        push!(st.polygons, collect(merged_dict[p]))
    end

    st.last_completed_picard_index = max_picard_index

    return st
end


@doc raw"""
    struct HDFPicardIndexStoragePreferences{T <: Integer}   

A struct holding preferences for Springer's classification using the HDF5 file format. It has four fields:

- `swmr :: Bool`: Whether to use single-reader-multiple-writer mode for HDF5.
    Defaults to `true`.
- `step_size :: Int`: The step size for multithreaded classification in
    terms of the picard index. Defaults to `10^4`.
- `maximum_picard_index :: Int`: The maximum picard index to be
    classified. Defaults to `10^7`.

"""
struct HDFPicardIndexStoragePreferences{T <: Integer}
    swmr :: Bool
    step_size :: Int
    maximum_picard_index :: Int
    
    HDFPicardIndexStoragePreferences{T}(
        swmr :: Bool = true,
        step_size :: Int = 10^4,
        maximum_picard_index :: Int = 10^7) where {T <: Integer} =
    new{T}(swmr, step_size, maximum_picard_index)

end


@doc raw"""
    mutable struct HDFPicardIndexStorage{T <: Integer} <: PicardIndexStorage{T}

A struct for managing classification results of Springer's classification of
lattice triangles using the HDF5 file format. It has the following fields:

- `preferences :: HDFPicardIndexStoragePreferences{T}`
- `file_path :: String`: The path of the HDF file to be generated.
- `last_completed_picard_index :: Int`: The last completed step of the classification. Initially, this will be `0`.
- `total_count :: Int`

"""
mutable struct HDFPicardIndexStorage{T <: Integer} <: PicardIndexStorage{T}
    preferences :: HDFPicardIndexStoragePreferences{T}
    file_path :: String
    last_completed_picard_index :: Int
    total_count :: Int

    function HDFPicardIndexStorage{T}(
        preferences :: HDFPicardIndexStoragePreferences{T},
        file_path :: String) where {T <: Integer}

        f = h5open(file_path, "cw"; swmr = preferences.swmr)
        f["numbers_of_polygons"] = zeros(Int, preferences.maximum_picard_index) 
        close(f)

        new{T}(preferences, file_path, 0, 0)

    end

    HDFPicardIndexStorage{T}(file_path :: String;
                         swmr :: Bool = true,
                         step_size :: Int = 10^4,
                         maximum_picard_index :: Int = 10^7) where {T <: Integer} =
    HDFPicardIndexStorage{T}(HDFPicardIndexStoragePreferences{T}(swmr, step_size, maximum_picard_index), file_path)

end


@doc raw"""
    classify_lattice_triangles_by_picard_index(st :: HDFPicardIndexStorage{T}, max_picard_index :: T; logging :: Bool = false) where {T <: Integer}

Perform Springer's classification of lattice triangles up go
`max_picard_index`, storing the results in an HDF5 file

"""
function classify_lattice_triangles_by_picard_index(st :: HDFPicardIndexStorage{T}, max_picard_index :: T; logging :: Bool = false) where {T <: Integer}


    step_size = st.preferences.step_size
    p0 = st.last_completed_picard_index
    N = max_picard_index - p0
    number_of_steps = N ÷ step_size + 1

    for b = 1 : number_of_steps

        p_min, p_max = p0 + (b-1)*step_size + 1, p0 + min(b*step_size, N)
        logging && @info "[p = $p_min : $p_max]: Beginning classification"

        dicts = Dict{T, Set{RationalPolygon{T,3,6}}}[]
        for t = 1 : Threads.nthreads()
            push!(dicts, Dict{T, Set{RationalPolygon{T,3,6}}}())
        end

        # Perform the classification
        Threads.@threads for p = p_min : p_max
            t = Threads.threadid()
            dicts[t][p] = classify_lattice_triangles_by_picard_index(p)
        end
        merged_dict = merge!(dicts...)

        new_count = sum([length(v) for (k,v) in merged_dict])
        logging && @info "[p = $p_min : $p_max]: Classfication complete. Found $new_count"

        # Write triangles to the HDF5 file
        f = h5open(st.file_path, "r+"; swmr = st.preferences.swmr)
        for p = p_min : p_max
            write_polygon_dataset(f, "pi$p", collect(merged_dict[p]))
            f["numbers_of_polygons"][p] = length(merged_dict[p])
            st.total_count += length(merged_dict[p])
        end
        close(f)

        logging && @info "[p = $p_min : $p_max]: Writeout complete. Running total: $(st.total_count)"

        st.last_completed_picard_index = p_max

    end

    return st
end
