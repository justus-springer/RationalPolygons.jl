@doc raw"""
    degree(w :: SVector{3})

Return the degree of a fake weighted projective plane with given fake
weight vector `w`. This is given by the formula `sum(w)^2 // prod(w)`,
see for instance Proposition 3.7 of [HaKi24](@cite).

"""
degree(w :: SVector{3}) = sum(w)^2 // prod(w)


@doc raw"""
    fake_weight_vector(P :: RationalPolygon{T,3}) where {T <: Integer}

Return the fake weight vector of the fake weighted projective plane
associated to the LDP triangle `P`.

"""
fake_weight_vector(P :: RationalPolygon{T,3}) where {T <: Integer} =
SVector{3,T}(multiplicity(P,2), multiplicity(P,3), multiplicity(P,1))


@doc raw"""
    n_step_mutations(us :: Set{SVector{3,T}}, depth :: Int) where {T <: Integer}   

Given solution triples `us` of a squared Markov type equation, return all
solution triples that can be obtained from `us` by applying at most `depth`
many mutations. This function returns a vector of length `depth+1`, containing
for each level `i = 0 : depth` the set of solution triples after exactly `i`
mutation steps.

# Example

Get all solutions to the squared Markov type equation ``9xyz=(x+y+z)^2`` by
starting with the initial solution ``(1,1,1)`` and applying at most four
mutations. Compare with ``T(9)`` of Remark 2.5 of [HaKi24](@cite).

```jldoctest
julia> us = Set([SVector{3,Int}(1,1,1)])
Set{SVector{3, Int64}} with 1 element:
  [1, 1, 1]

julia> n_step_mutations(us, 4)
5-element Vector{Set{SVector{3, Int64}}}:
 Set([[1, 1, 1]])
 Set([[1, 1, 4]])
 Set([[1, 4, 25]])
 Set([[4, 25, 841], [1, 25, 169]])
 Set([[25, 169, 37636], [25, 841, 187489], [1, 169, 1156], [4, 841, 28561]])
```

"""
function n_step_mutations(us :: Set{SVector{3,T}}, depth :: Int) where {T <: Integer}
    res = Set{SVector{3,T}}[]
    push!(res, us)
    for i = 1 : depth
        next_triples = Set{SVector{3,T}}()
        for u ∈ last(res)
            push!(next_triples, SVector{3,T}(u[2], u[3], (u[2]+u[3])^2 ÷ u[1]))
            push!(next_triples, SVector{3,T}(u[1], u[3], (u[1]+u[3])^2 ÷ u[2]))
        end
        push!(res, next_triples)
    end
    return res
end


@doc raw"""
    initial_triple(:: Val{A}, T :: Type{<:Integer} = BigInt)

Return the unique initial triple of the squared Markov type equation ``Axyz =
(x+y+z)^2``. Allowed values of `A` are 9, 8, 6 and 5.

# Example

The unique initial triples for all allowed values of `A`, see Theorem 2.2 of
[HaKi24](@cite).

```jldoctest
julia> [initial_triple(Val(A)) for A in [9,8,6,5]]
4-element Vector{SVector{3, BigInt}}:
 [1, 1, 1]
 [1, 1, 2]
 [1, 2, 3]
 [1, 4, 5]
```

"""
initial_triple(:: Val{9}, T :: Type{<:Integer} = BigInt) = SVector{3,T}(1,1,1)
initial_triple(:: Val{8}, T :: Type{<:Integer} = BigInt) = SVector{3,T}(1,1,2)
initial_triple(:: Val{6}, T :: Type{<:Integer} = BigInt) = SVector{3,T}(1,2,3)
initial_triple(:: Val{5}, T :: Type{<:Integer} = BigInt) = SVector{3,T}(1,4,5)


@doc raw"""
    classify_squared_markov_type_equation_solutions(::Val{A}, depth :: Int, T :: Type{<:Integer} = BigInt) where {A}

Return all solutions to the squared Markov type equation ``Axyz = (x+y+z)^2``
by starting with the initial solution and applying at most `depth` many
mutations. Allowed values of `A` are 9, 8, 6 and 5.

# Example

All solution triples to ``8xyz = (x+y+z)^2`` with depth at most three.

```jldoctest
julia> classify_squared_markov_type_equation_solutions(Val(8),3)
4-element Vector{Set{SVector{3, BigInt}}}:
 Set([[1, 1, 2]])
 Set([[1, 2, 9]])
 Set([[2, 9, 121], [1, 9, 50]])
 Set([[9, 50, 3481], [2, 121, 1681], [1, 50, 289], [9, 121, 8450]])
```

"""
function classify_squared_markov_type_equation_solutions(::Val{A}, depth :: Int, T :: Type{<:Integer} = BigInt) where {A}
    us = Set{SVector{3,T}}()
    push!(us, initial_triple(Val(A), T))
    return n_step_mutations(us, depth)
end


@doc raw"""
    adjust_triple(::Val{A}, u :: SVector{3,T}) where {T <: Integer}
    adjust_triple(u :: SVector{3,T}) where {T <: Integer}

Reorder a solution triple of a squared Markov type equation to make it
adjusted, according to Definition 3.14 of [HaKi24](@cite). Allowed values of
`A` are 9, 8, 6 and 5. If `A` is not given, it is determined by calculating the
degree of `u`.

# Example

The solution triple ``(50,9,3481)`` for ``A=8`` is not adjusted, since the even
entry is not last.

```jldoctest
julia> adjust_triple(Val(8), SVector{3,Int}(50,9,3481))
3-element SVector{3, Int64} with indices SOneTo(3):
    9
 3481
   50
```

"""
function adjust_triple(::Val{9}, u :: SVector{3,T}) where {T <: Integer}
    return sort(u)
end

function adjust_triple(::Val{8}, u :: SVector{3,T}) where {T <: Integer}
    i = findfirst(x -> x % 2 == 0, u)
    j, k = filter(x -> x ≠ i, 1 : 3)
    if u[j] ≤ u[k]
        return SVector{3,T}(u[j],u[k],u[i])
    else
        return SVector{3,T}(u[k],u[j],u[i])
    end
end

function adjust_triple(::Val{6}, u :: SVector{3,T}) where {T <: Integer}
    i = findfirst(x -> x % 2 == 0, u)
    j = findfirst(x -> x % 3 == 0, u)
    k = first(filter(x -> x ≠ i && x ≠ j, 1 : 3))
    return SVector{3,T}(u[k],u[i],u[j])
end

function adjust_triple(::Val{5}, u :: SVector{3,T}) where {T <: Integer}
    i = findfirst(x -> x % 5 == 0, u)
    j, k = filter(x -> x ≠ i, 1 : 3)
    if u[j] ≤ u[k]
        return SVector{3,T}(u[j],u[k],u[i])
    else
        return SVector{3,T}(u[k],u[j],u[i])
    end
end

function adjust_triple(u :: SVector{3,T}) where {T <: Integer}
    K = degree(u)
    denominator(K) == 1 || error("degree is not integral")
    return adjust_triple(Val(Int(numerator(K))), u)
end


@doc raw"""
    fake_weight_vectors_to_triangles(us :: Set{SVector{3,T}}, μ :: T) where {T <: Integer}

Given a set of triples `us` sharing the same integral degree `A` and an integer
`μ`, return all LDP triangles having fake weight vector `μ * u`. Allowed values
of `A` are 9, 8, 6 and 5. Moreover, `μ` must be a divisor of `A`. The resulting
triangles will have degree `A ÷ μ` and multiplicity `μ`.

"""
function fake_weight_vectors_to_triangles(us :: Set{SVector{3,T}}, μ :: T) where {T <: Integer}
    res = Set{RationalPolygon{T,3}}()
    for w ∈ us
        w0,w1,w2 = adjust_triple(w)
        _, x, y = gcdx(w1, w2)
        for k = 1 : μ
            a = -w0*(x - k*w2)
            b = -w0*(y + k*w1)
            P = RationalPolygon(SMatrix{2,3,T}(1,0,a,μ*w2,b,-μ*w1), one(T))
            is_primitive(P) || continue
            push!(res, unimodular_normal_form(P))
        end
    end
    return res
end


@doc raw"""
    classify_lattice_triangles_integral_degree(::Val{K}, ::Val{μ}, depth :: Int, T :: Type{<:Integer} = BigInt) where {K, μ}

Return all LDP triangles (= fake weighted projective planes) with integral
degree `K` and class group torsion order `μ`, up to a given depth in the Markov
tree. Essentially, this returns the triangles of the series (K-μ-*) according
to the notation of Theorem 1.1 of [HaKi24](@cite). Allowed values of `K` are
1, 2, 3, 4, 5, 6, 8 and 9. Allowed values of `μ` are all integers such that
`μ * K` is 5, 6, 8 or 9.

# Example

Compute all LDP triangles of degree 1 and class group torsion order 9, up to
Markov depth 5. These consist of the three series (1-9-2), (1-9-5) and (1-9-8)
from Theorem 1.1 of [HaKi24](@cite). Note that for depths zero and one (which
correspond to solution triples (1,1,1) and (1,1,4) in the Markov tree), the
series overlap, hence there are only one resp. two triangles in this case,
which is also mentioned in the Theorem. After that, the number of triangles
doubles with each additional step.

```jldoctest
julia> Pss = classify_lattice_triangles_integral_degree(Val(1), Val(9), 5);

julia> length.(Pss)
6-element Vector{Int64}:
  1
  2
  3
  6
 12
 24
```

"""
function classify_lattice_triangles_integral_degree(::Val{K}, ::Val{μ}, depth :: Int, T :: Type{<:Integer} = BigInt) where {K, μ}
    A = K * μ
    wss = classify_squared_markov_type_equation_solutions(Val(A), depth, T)
    return [fake_weight_vectors_to_triangles(ws, T(μ)) for ws ∈ wss]
end


@doc raw"""
    classify_lattice_triangles_integral_degree(::Val{K}, depth :: Int, T :: Type{<:Integer} = BigInt) where {K}

Return all LDP triangles (= fake weighted projective planes) with integral
degree `K`, up to a given depth in the Markov tree. Essentially, this returns
the triangles of the series (K - * - *) according to the notation of Theorem 1.1 of
[HaKi24](@cite). Allowed values of `K` are 1, 2, 3, 4, 5, 6, 8 and 9.

# Example

Print the numbers of LDP triangles with integral degree `K`, for all possible
values of `K`, up to depth 10.

```jldoctest
julia> [length.(classify_lattice_triangles_integral_degree(Val(K), 10)) for K in [1,2,3,4,5,6,8,9]]
8-element Vector{Vector{Int64}}:
 [10, 18, 35, 70, 140, 280, 560, 1120, 2240, 4480, 8960]
 [4, 6, 12, 24, 48, 96, 192, 384, 768, 1536, 3072]
 [2, 3, 5, 10, 20, 40, 80, 160, 320, 640, 1280]
 [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
 [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
 [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
 [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
 [1, 1, 1, 2, 4, 8, 16, 32, 64, 128, 256]
```

"""
classify_lattice_triangles_integral_degree(::Val{9}, depth :: Int, T :: Type{<:Integer} = BigInt) =
classify_lattice_triangles_integral_degree(Val(9), Val(1), depth, T)

classify_lattice_triangles_integral_degree(::Val{8}, depth :: Int, T :: Type{<:Integer} = BigInt) =
classify_lattice_triangles_integral_degree(Val(8), Val(1), depth, T)

classify_lattice_triangles_integral_degree(::Val{6}, depth :: Int, T :: Type{<:Integer} = BigInt) =
classify_lattice_triangles_integral_degree(Val(6), Val(1), depth, T)

classify_lattice_triangles_integral_degree(::Val{5}, depth :: Int, T :: Type{<:Integer} = BigInt) =
classify_lattice_triangles_integral_degree(Val(5), Val(1), depth, T)

classify_lattice_triangles_integral_degree(::Val{4}, depth :: Int, T :: Type{<:Integer} = BigInt) =
classify_lattice_triangles_integral_degree(Val(4), Val(2), depth, T)

function classify_lattice_triangles_integral_degree(::Val{3}, depth :: Int, T :: Type{<:Integer} = BigInt)
    Ps2 = classify_lattice_triangles_integral_degree(Val(3), Val(2), depth, T)
    Ps3 = classify_lattice_triangles_integral_degree(Val(3), Val(3), depth, T)
    for i = 1 : length(Ps2)
        union!(Ps2[i], Ps3[i])
    end
    return Ps2
end

function classify_lattice_triangles_integral_degree(::Val{2}, depth :: Int, T :: Type{<:Integer} = BigInt)
    Ps3 = classify_lattice_triangles_integral_degree(Val(2), Val(3), depth, T)
    Ps4 = classify_lattice_triangles_integral_degree(Val(2), Val(4), depth, T)
    for i = 1 : length(Ps3)
        union!(Ps3[i], Ps4[i])
    end
    return Ps3
end

function classify_lattice_triangles_integral_degree(::Val{1}, depth :: Int, T :: Type{<:Integer} = BigInt)
    Ps5 = classify_lattice_triangles_integral_degree(Val(1), Val(5), depth, T)
    Ps6 = classify_lattice_triangles_integral_degree(Val(1), Val(6), depth, T)
    Ps8 = classify_lattice_triangles_integral_degree(Val(1), Val(8), depth, T)
    Ps9 = classify_lattice_triangles_integral_degree(Val(1), Val(9), depth, T)
    for i = 1 : length(Ps5)
        union!(Ps5[i], Ps6[i], Ps8[i], Ps9[i])
    end
    return Ps5
end


@doc raw"""
    classify_lattice_triangles_integral_degree(depth :: Int, T :: Type{<:Integer} = BigInt)

Return all LDP triangles (= fake weighted projective planes) with integral
degree, up to a given depth in the Markov tree. This returns the union of all
24 series classified in Theorem 1.1 of [HaKi24](@cite).

"""
function classify_lattice_triangles_integral_degree(depth :: Int, T :: Type{<:Integer} = BigInt)
    Pss = [classify_lattice_triangles_integral_degree(Val(K), depth, T) for K in [1,2,3,4,5,6,8,9]]
    return [union!([Pss[K][i] for K = 1 : 8]...) for i = 1 : depth+1]
end
