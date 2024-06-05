
floor_k_rational(k :: T, x :: Real) where {T <: Integer} =
floor(T, k * x) // k

ceil_k_rational(k :: T, x :: Real) where {T <: Integer} =
ceil(T, k * x) // k

