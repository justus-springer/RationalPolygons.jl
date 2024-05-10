
function reset_attributes!(x)
    try
        empty!(x.__attrs)
    finally
        return x
    end
end

reset_attributes!(xs :: AbstractVector) = map(reset_attributes!, xs)


floor_k_rational(k :: T, x :: Real) where {T <: Integer} =
floor(T, k * x) // k

ceil_k_rational(k :: T, x :: Real) where {T <: Integer} =
ceil(T, k * x) // k

