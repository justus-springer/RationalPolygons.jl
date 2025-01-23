export solve_unique_integer_solution

function solve_unique_integer_solution(A :: SMatrix{N,N,T}, b :: SVector{N,T}) where {N, T <: Integer}
    H,U = hnfc(A)
    x = T[]
    for i = 1 : N
        d = b[i] - sum(T[H[i,k]*x[k] for k = 1 : i-1])
        d % H[i,i] == 0 || error("No integral solution")
        push!(x, d ÷ H[i,i])
    end
    return SVector{N,T}(U*x)
end
