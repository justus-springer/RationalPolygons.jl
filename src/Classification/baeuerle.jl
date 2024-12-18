@doc raw"""
    unit_fraction_partitions_length_three(ι :: T) where {T <: Integer}

Return all triples `(a,b,c)` such that `1//ι = 1//a + 1//b + 1//c` Andreas`a ≤ b ≤ c`. See also A004194 on oeis.

# Example

```jldoctest
julia> unit_fraction_partitions_length_three(2)
10-element Vector{Tuple{Int64, Int64, Int64}}:
 (3, 7, 42)
 (3, 8, 24)
 (3, 9, 18)
 (3, 10, 15)
 (3, 12, 12)
 (4, 5, 20)
 (4, 6, 12)
 (4, 8, 8)
 (5, 5, 10)
 (6, 6, 6)
```

"""
function unit_fraction_partitions_length_three(ι :: T) where {T <: Integer}
    res = Tuple{T,T,T}[]
    for a = ι + 1 : 3ι
        for b = max(fld(a*ι, a - ι) + 1, a) : fld(2*a*ι, a - ι)
            (c,r) = divrem(ι*a*b, a*b - a*ι - b*ι)
            r == 0 || continue
            push!(res, (a,b,c))
        end
    end
    return res
end


@doc raw"""
    classify_lattice_triangles_by_gorenstein_index(ι :: T) where {T <: Integer}

Return all lattice triangles with gorenstein index `ι`.

# Example

There are five lattice triangles with gorenstein index one:

```jldoctest
julia> classify_lattice_triangles_by_gorenstein_index(1)
Set{RationalPolygon{Int64, 3}} with 5 elements:
  Rational polygon of rationality 1 with 3 vertices.
  Rational polygon of rationality 1 with 3 vertices.
  Rational polygon of rationality 1 with 3 vertices.
  Rational polygon of rationality 1 with 3 vertices.
  Rational polygon of rationality 1 with 3 vertices.
```

"""
function classify_lattice_triangles_by_gorenstein_index(ι :: T) where {T <: Integer}
    res = Set{RationalPolygon{T,3}}()

    for (a0, a1, a2) in unit_fraction_partitions_length_three(ι)

        # Compute the weight vector from the unit fraction partition
        g = gcd(a1 * a2, a0 * a2, a0 * a1)
        w0, w1, w2 = a1 * a2 ÷ g, a0 * a2 ÷ g, a0 * a1 ÷ g

        # Check if it's well formed
        (gcd(w0, w1) == 1 && gcd(w0, w2) == 1 && gcd(w1, w2) == 1) || continue

        γ_max = gcd(a0,a1)
        for γ = 1 : γ_max
            # `γ` must be a divisor of `gcd(a0,a1)`
            γ_max % γ == 0 || continue
            for α = 0 : γ-1
                gcd(α, γ) == 1 || continue
                (w0 + w1 * α) % w2 == 0 || continue
                (w1 * γ) % w2 == 0 || continue
                β = -(w0 + w1 * α) ÷ w2
                δ = - w1 * γ ÷ w2
                gcd(β, δ) == 1 || continue

                P = RationalPolygon(SMatrix{2,3,T,6}(1,0,α,γ,β,δ), 1)
                gorenstein_index(P) == ι || continue
                push!(res, unimodular_normal_form(P))
            end
        end
    end
    return res
end


@doc raw"""
    abstract type BaeuerleStorage{T <: Integer} end

Abstract supertype of `InMemoryBaeuerleStorage` and `HDFBaeuerleStorage`.

"""
abstract type BaeuerleStorage{T <: Integer} end


@doc raw"""
    mutable struct InMemoryBaeuerleStorage{T <: Integer} <: BaeuerleStorage{T}

A struct holding classification results of Baeuerle's classification of lattice
triangles by gorenstein index.

"""
mutable struct InMemoryBaeuerleStorage{T <: Integer} <: BaeuerleStorage{T}
    polygons :: Vector{Vector{RationalPolygon{T,3,6}}}
    last_completed_gorenstein_index :: T

    InMemoryBaeuerleStorage{T}() where {T <: Integer} = new{T}(Vector{RationalPolygon{T,3,6}}[], 0)

end


@doc raw"""
    classify_lattice_triangles_by_gorenstein_index(st :: InMemoryBaeuerleStorage{T}, max_gorenstein_index :: T) where {T <: Integer}

Perform Bäuerle's classification of lattice triangles up go
`max_gorenstein_index`, storing the results in memory.

# Example

Reproduce Bäuerle's original classification up to gorenstein index 1000, see
Theorem 1.4 of [Bae23](@cite).

```jldoctest
julia> st = InMemoryBaeuerleStorage{Int}()
InMemoryBaeuerleStorage{Int64}(Vector{RationalPolygon{Int64, 3, 6}}[], 0)

julia> classify_lattice_triangles_by_gorenstein_index(st, 1000);

julia> sum(length.(st.polygons))
2992229
```

"""
function classify_lattice_triangles_by_gorenstein_index(st :: InMemoryBaeuerleStorage{T}, max_gorenstein_index :: T) where {T <: Integer}
    dicts = Dict{T, Set{RationalPolygon{T,3,6}}}[]
    for t = 1 : Threads.nthreads()
        push!(dicts, Dict{T, Set{RationalPolygon{T,3,6}}}())
    end

    ι0 = st.last_completed_gorenstein_index + 1
    Threads.@threads for ι = ι0 : max_gorenstein_index
        t = Threads.threadid()
        dicts[t][ι] = classify_lattice_triangles_by_gorenstein_index(ι)
    end

    merged_dict = merge!(dicts...)
    for ι = ι0 : max_gorenstein_index
        push!(st.polygons, collect(merged_dict[ι]))
    end

    st.last_completed_gorenstein_index = max_gorenstein_index

    return st
end


@doc raw"""
    struct HDFBaeuerleStoragePreferences{T <: Integer}   

A struct holding preferences for Baeuerle's classification using the HDF5 file format. It has four fields:

- `swmr :: Bool`: Whether to use single-reader-multiple-writer mode for HDF5.
    Defaults to `true`.
- `step_size :: Int`: The step size for multithreaded classification in
    terms of the gorenstein index. Defaults to `100`.
- `maximum_gorenstein_index :: Int`: The maximum gorenstein index to be
    classified. Defaults to `10^5`.

"""
struct HDFBaeuerleStoragePreferences{T <: Integer}
    swmr :: Bool
    step_size :: Int
    maximum_gorenstein_index :: Int
    
    HDFBaeuerleStoragePreferences{T}(
        swmr :: Bool = true,
        step_size :: Int = 100,
        maximum_gorenstein_index :: Int = 10^5) where {T <: Integer} =
    new{T}(swmr, step_size, maximum_gorenstein_index)

end


@doc raw"""
    mutable struct HDFBaeuerleStorage{T <: Integer} <: BaeuerleStorage{T}

A struct for managing classification results of Baeuerle's classification of
lattice triangles using the HDF5 file format. It has the following fields:

- `preferences :: HDFBaeuerleStoragePreferences{T}`
- `file_path :: String`: The path of the HDF file to be generated.
- `last_completed_gorenstein_index :: Int`: The last completed step of the classification. Initially, this will be `3`.
- `total_count :: Int`

"""
mutable struct HDFBaeuerleStorage{T <: Integer} <: BaeuerleStorage{T}
    preferences :: HDFBaeuerleStoragePreferences{T}
    file_path :: String
    last_completed_gorenstein_index :: Int
    total_count :: Int

    function HDFBaeuerleStorage{T}(
        preferences :: HDFBaeuerleStoragePreferences{T},
        file_path :: String) where {T <: Integer}

        f = h5open(file_path, "cw"; swmr = preferences.swmr)
        f["numbers_of_polygons"] = zeros(Int, preferences.maximum_gorenstein_index) 

        close(f)

        new{T}(preferences, file_path, 0, 0)

    end

    HDFBaeuerleStorage{T}(file_path :: String;
                         swmr :: Bool = true,
                         step_size :: Int = 100,
                         maximum_gorenstein_index :: Int = 10^5) where {T <: Integer} =
    HDFBaeuerleStorage{T}(HDFBaeuerleStoragePreferences{T}(swmr, step_size, maximum_gorenstein_index), file_path)

end


@doc raw"""
    classify_lattice_triangles_by_gorenstein_index(st :: HDFBaeuerleStorage{T}, max_gorenstein_index :: T; logging :: Bool = false) where {T <: Integer}

Perform Bäuerle's classification of lattice triangles up go
`max_gorenstein_index`, storing the results in an HDF5 file

"""
function classify_lattice_triangles_by_gorenstein_index(st :: HDFBaeuerleStorage{T}, max_gorenstein_index :: T; logging :: Bool = false) where {T <: Integer}


    step_size = st.preferences.step_size
    ι0 = st.last_completed_gorenstein_index
    N = max_gorenstein_index - ι0
    number_of_steps = N ÷ step_size + 1

    for b = 1 : number_of_steps

        ι_min, ι_max = ι0 + (b-1)*step_size + 1, ι0 + min(b*step_size, N)
        logging && @info "[ι = $ι_min : $ι_max]: Beginning classification"

        dicts = Dict{T, Set{RationalPolygon{T,3,6}}}[]
        for t = 1 : Threads.nthreads()
            push!(dicts, Dict{T, Set{RationalPolygon{T,3,6}}}())
        end

        # Perform the classification
        Threads.@threads for ι = ι_min : ι_max
            t = Threads.threadid()
            dicts[t][ι] = classify_lattice_triangles_by_gorenstein_index(ι)
        end
        merged_dict = merge!(dicts...)

        new_count = sum([length(v) for (k,v) in merged_dict])
        logging && @info "[ι = $ι_min : $ι_max]: Classfication complete. Found $new_count"

        # Write triangles to the HDF5 file
        f = h5open(st.file_path, "r+"; swmr = st.preferences.swmr)
        for ι = ι_min : ι_max
            write_polygon_dataset(f, "gi$ι", collect(merged_dict[ι]))
            f["numbers_of_polygons"][ι] = length(merged_dict[ι])
            st.total_count += length(merged_dict[ι])
        end
        close(f)

        logging && @info "[ι = $ι_min : $ι_max]: Writeout complete. Running total: $(st.total_count)"

        st.last_completed_gorenstein_index = ι_max

    end

    return st
end
