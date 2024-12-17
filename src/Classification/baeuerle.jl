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
