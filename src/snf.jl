#=

This file provides an elementary implementation of the Smith normal form for
integer matrices. The implementation is heavily inspired by the one used in the
NormalForms.jl package written by Brandon Flores, which itself uses code from
the HermiteNormalForms.jl package authored by Yingbo Ma and Chris Elrod. These
repositories can be found at:

https://github.com/brainandforce/NormalForms.jl
https://github.com/YingboMa/HermiteNormalForm.jl

=#
  
@doc raw"""
    modify_rows!(
        A :: AbstractMatrix{T},
        i1 :: Integer,
        i2 :: Integer,
        α :: T, β :: T, γ :: T, δ :: T) where {T <: Integer}

A modification of two rows that occurs in the Smith normal form algorithm. It
replaces the `i1`-th row with `α` times the original `i1`-th row plus `β` times
the original `i2`-th row, and the `i2`-th row with `γ` times the original
`i1`-th row plus `δ` times the original `i2`-th row. This is a unimodular
transformation if and only if the matrix `[α β ; γ δ]` is unimodular, i.e.
`α * δ - β * γ` is plus or minus one.

"""
function modify_rows!(
        A :: AbstractMatrix{T},
        i1 :: Integer,
        i2 :: Integer,
        α :: T, β :: T, γ :: T, δ :: T) where {T <: Integer}

    for j in axes(A, 2)
        x, y = A[i1,j], A[i2,j]
        A[i1,j] = α * x + β * y
        A[i2,j] = γ * x + δ * y

    end

end


@doc raw"""
    modify_columns!(
        A :: AbstractMatrix{T},
        j1 :: Integer,
        j2 :: Integer,
        α :: T, β :: T, γ :: T, δ :: T) where {T <: Integer}

A modification of two columns that occurs in the Smith normal form algorithm.
It replaces the `j1`-th column with `α` times the original `j1`-th column plus
`β` times the original `j2`-th column, and the `j2`-th column with `γ` times
the original `j1`-th column plus `δ` times the original `j2`-th column. This is
a unimodular transformation if and only if the matrix `[α β ; γ δ]` is
unimodular, i.e. `α * δ - β * γ` is plus or minus one.

"""
function modify_columns!(
        A :: AbstractMatrix{T},
        j1 :: Integer,
        j2 :: Integer,
        α :: T, β :: T, γ :: T, δ :: T) where {T <: Integer}

    for i in axes(A, 1)
        x, y = A[i,j1], A[i,j2]
        A[i,j1] = α * x + β * y
        A[i,j2] = γ * x + δ * y
    end

end


@doc raw"""
    add_multiple_to_row!(
        A :: AbstractMatrix{T},
        i1 :: Integer,
        i2 :: Integer,
        t :: T) where {T <: Integer}

Add the `t`-fold of the `i1`-th row to the `i2`-th row.

"""
function add_multiple_to_row!(
        A :: AbstractMatrix{T},
        i1 :: Integer,
        i2 :: Integer,
        t :: T) where {T <: Integer}

    for j in axes(A, 2)
        A[i2, j] += t * A[i1, j]
    end

end


@doc raw"""
    add_multiple_to_column!(
        A :: AbstractMatrix{T},
        j1 :: Integer,
        j2 :: Integer,
        t :: T) where {T <: Integer}

Add the `t`-fold of the `j1`-th column to the `j2`-th column.

"""
function add_multiple_to_column!(
        A :: AbstractMatrix{T},
        j1 :: Integer,
        j2 :: Integer,
        t :: T) where {T <: Integer}

    for i in axes(A, 1)
        A[i, j2] += t * A[i, j1]
    end

end


@doc raw"""
    invert_row!(
        A :: AbstractMatrix{T},
        i :: Integer) where {T <: Integer}

Multiply the `i`-th row of `A` by minus one.

"""
function invert_row!(
        A :: AbstractMatrix{T},
        i :: Integer) where {T <: Integer}

    for j in axes(A,2)
        A[i,j] *= -1
    end

end


@doc raw"""
    invert_column!(
        A :: AbstractMatrix{T},
        j :: Integer) where {T <: Integer}

Multiply the `j`-th column of `A` by minus one.

"""
function invert_column!(
        A :: AbstractMatrix{T},
        j :: Integer) where {T <: Integer}

    for i in axes(A,1)
        A[i,j] *= -1
    end

end


@doc raw"""
    zero_row_after!(
        A :: AbstractMatrix{T},
        V :: AbstractMatrix{T},
        i :: Integer,
        j :: Integer) where {T <: Integer}

Use unimodular transformations to achieve `A[i,k] == 0` for all j+1 ≤ k ≤ n,
where `n` is the number of columns of `A`.

"""
function zero_row_after!(
        A :: AbstractMatrix{T},
        V :: AbstractMatrix{T},
        i :: Integer,
        j :: Integer) where {T <: Integer}

    # First, we go through all j+1 ≤ k ≤ n such that `A[i,j]` (the pivot) does
    # not divide `A[i,k]`. We then apply unimodular transformations to replace
    # `A[i,j]` with the gcd of `A[i,j]` and `A[i,k]`.

    for k = j+1 : maximum(axes(A,2))
        if A[i,k] % A[i,j] != 0
            d, α, β = gcdx(A[i,j], A[i,k])
            γ, δ = A[i,j] ÷ d, A[i,k] ÷ d
            modify_columns!(A, j, k, α, β, -δ, γ)
            modify_columns!(V, j, k, α, β, -δ, γ)
        end
    end

    # Now, we have guaranteed `A[i,k] % A[i,j] == 0` for all j+1 ≤ k ≤ n.
    # Hence we can multiples of the `j`-th column to the `k`-th column
    # to achieve `A[i,k] == 0` for all j+1 ≤ k ≤ n.
    
    for k = j+1 : maximum(axes(A,2))
        x = A[i,k] ÷ A[i,j]
        add_multiple_to_column!(A, j, k, -x)
        add_multiple_to_column!(V, j, k, -x)
    end

end


@doc raw"""
    zero_column_after!(
        A :: AbstractMatrix{T},
        U :: AbstractMatrix{T},
        i :: Integer,
        j :: Integer) where {T <: Integer}

Use unimodular transformations to achieve `A[k,j] == 0` for all i+1 ≤ k ≤ m,
where `m` the number of rows of `A`.

"""
function zero_column_after!(
        A :: AbstractMatrix{T},
        U :: AbstractMatrix{T},
        i :: Integer,
        j :: Integer) where {T <: Integer}

    # First, we go through all i+1 ≤ k ≤ m such that `A[i,j]` (the pivot) does
    # not divide `A[k,j]`. We then apply unimodular transformations to replace
    # `A[i,j]` with the gcd of `A[i,j]` and `A[k,j]`.

    for k = i+1 : maximum(axes(A,1))
        if A[k,j] % A[i,j] != 0
            d, α, β = gcdx(A[i,j], A[k,j])
            γ, δ = A[i,j] ÷ d, A[k,j] ÷ d
            modify_rows!(A, i, k, α, β, -δ, γ)
            modify_rows!(U, i, k, α, β, -δ, γ)
        end
    end

    # Now we have `A[k,j] % A[i,j] == 0` for all i+1 ≤ k ≤ m.
    # Hence we can multiples of the `i`-th row to the `k`-th row
    # to achieve `A[k,j] == 0` for all i+1 ≤ k ≤ m.
    
    for k = i+1 : maximum(axes(A,1))
        x = A[k,j] ÷ A[i,j]
        add_multiple_to_row!(A, i, k, -x)
        add_multiple_to_row!(U, i, k, -x)
    end

end


is_zero_row(A :: AbstractMatrix, i :: Integer) = all(iszero(A[i,:]))
is_zero_row_after(A::AbstractMatrix, i::Integer, j ::Integer) = all(iszero(A[i,j+1:end]))


is_zero_column(A :: AbstractMatrix, j :: Integer) = all(iszero(A[:,j]))
is_zero_column_after(A::AbstractMatrix, i::Integer, j ::Integer) = all(iszero(A[i+1:end,j]))


@doc raw"""
    zero_row_and_column_after!(
        A :: AbstractMatrix{T},
        U :: AbstractMatrix{T},
        V :: AbstractMatrix{T},
        i :: Integer,
        j :: Integer) where {T <: Integer}

Apply unimodular transformations to `A` to achieve `A[i,k] == 0` for all
`j+1 ≤ k ≤ n` and `A[k,j] == 0` for all `i+1 ≤ k ≤ m`, where `m` and `n` are
the numbers of rows and columns of `A` respectively.

"""
function zero_row_and_column_after!(
        A :: AbstractMatrix{T},
        U :: AbstractMatrix{T},
        V :: AbstractMatrix{T},
        i :: Integer,
        j :: Integer) where {T <: Integer}

    while !is_zero_row_after(A, i, j) || !is_zero_column_after(A, i, j)
        zero_row_after!(A, V, i, j)
        zero_column_after!(A, U, i, j)
    end

end


function snf_initialize_transformation_matrices(A :: AbstractMatrix{T}) where {T <: Integer}
    m, n = size(A,1), size(A,2)
    U = typeof(A)(collect(UniformScaling(one(T))(m)))
    V = typeof(A)(collect(UniformScaling(one(T))(n)))
    return U, V
end

function snf_initialize_transformation_matrices(A :: StaticMatrix{M,N,T}) where {M, N, T <: Integer}
    U = MMatrix{M,M}(collect(UniformScaling(one(T))(M)))
    V = MMatrix{N,N}(collect(UniformScaling(one(T))(N)))
    return U, V
end


@doc raw"""
    snf_with_transform!(A :: AbstractMatrix{T}) where {T <: Integer}

Apply unimodular transformations to turn `A` into Smith normal form. This
function modifies its input. It returns a triple `(A, U, V)`, where `U` is the
transformation matrix for the rows and `V` is the transformation matrix of the
columns. They satisfy `U * A * V == S`, where `S` is the Smith normal form of
`A`.

"""
function snf_with_transform!(A :: AbstractMatrix{T}) where {T <: Integer}
    
    # Let m be the number of rows and n the number of columns.
    m, n = size(A,1), size(A,2)

    # Initialize the transformation matrices as identity matrices.
    U, V = snf_initialize_transformation_matrices(A)

    # We loop through the smaller dimension of `A`.
    for k in minimum(axes(A))

        # We would like to have a nonzero pivot at position (k,k). If the
        # `k`-th column is zero, we first have to swap the `k`-th column with
        # some nonzero column.
        if is_zero_column(A, k)

            j = findfirst(j -> !is_zero_column(A, j), k+1 : n)

            # If `j` is nothing, this means there are non nonzero columns left,
            # hence we are done and can break out of the loop immediately.
            isnothing(j) && break

            # Otherwise, we can swap columns.
            swapcols!(A, k, k + j)
            swapcols!(V, k, k + j)
        end

        # Now we know that the `k`-th column is nonzero. However, we still
        # can't guarantee that `A[k,k]` is nonzero. If it is zero, we have
        # to find a nonzero entry `A[i,k]` and swap rows.
        if iszero(A[k,k])

            i = findfirst(i -> !iszero(A[i,k]), k+1 : m)

            # Note that `i` is can't be nothing, as we already know this column
            # contains some nonzero entry.
            
            swaprows!(A, k, k + i)
            swaprows!(U, k, k + i)
        end

        # Now we have `A[k,k] != 0`. Now we perform the main step of the
        # algorithm to achieve zero entries to the right and below of the
        # pivot.
        zero_row_and_column_after!(A, U, V, k, k)

        # We want to have strictly positive entries in the diagonal
        if A[k,k] < 0
            invert_column!(A, k)
            invert_column!(V, k)
        end

        # This loop ensures the divisibility condition between the diagonal
        # entries. If we find that some diagonal entry does not divide the next
        # one, we have to reconcile this by means of unimodular
        # transformations.
        for l = 1 : k-1
            if A[k-l+1,k-l+1] % A[k-l,k-l] != 0
                add_multiple_to_row!(A, k-l+1, k-l, one(T))
                add_multiple_to_row!(U, k-l+1, k-l, one(T))
                zero_row_and_column_after!(A, U, V, k-l, k-l)
            end
        end
    end

    return A, U, V

end


@doc raw"""
    snf_with_transform(A :: AbstractMatrix)

Non-modifying version of `snf_with_transform!`.

"""
snf_with_transform(A :: AbstractMatrix) = snf_with_transform!(deepcopy(A))

function snf_with_transform(A :: SMatrix)
    S, U, V = snf_with_transform!(convert(MMatrix,A))
    return convert(SMatrix, S), convert(SMatrix, U), convert(SMatrix, V)
end


@doc raw"""
    snf!(A :: AbstractMatrix)

Turn `A` into Smith normal form. This function modifies its input.

"""
snf!(A :: AbstractMatrix) = snf_with_transform!(A)[1]


@doc raw"""
    snf(A :: AbstractMatrix)

Non-modifying version of `snf!`.

"""
snf(A :: AbstractMatrix) = snf_with_transform(A)[1]


@doc raw"""
    elementary_divisors(A :: AbstractMatrix{T}) where {T <: Integer}

Return the elementary_divisors divisors of `A`, i.e. the nonzero diagonal
entries in the Smith normal form.

"""
function elementary_divisors(A :: AbstractMatrix{T}) where {T <: Integer}
    res = T[]
    S = snf(A)
    for k = minimum(axes(A))
        if S[k,k] > 0
            push!(res, S[k,k])
        end
    end
    return res
end


@doc raw"""
    is_snf(A :: AbstractMatrix)

Check whether the matrix `A` is in Smith normal form.

"""
function is_snf(A :: AbstractMatrix)
    m, n = maximum(axes(A,1)), maximum(axes(A,2))
    r = length(filter(i -> A[i,i] ≠ 0, 1 : min(n,m)))

    for i in axes(A,1)
        for j in axes(A,2)
            if i == j && i <= r
                A[i,i] > 0 || return false
            else
                A[i,j] == 0 || return false
            end
        end
    end

    # check the divisibility condition
    all(k -> A[k,k] % A[k-1,k-1] == 0, 2 : r) || return false

    return true

end
