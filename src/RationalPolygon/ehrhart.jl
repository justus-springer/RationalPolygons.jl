@doc raw"""
    is_periodic(v :: Vector, k :: Int)

Check if a vector `v` is `k`-periodic, i.e. `v[i] == v[mod(i+k,1:n)]` for all
`i`, where `n = length(v)`.

"""
function is_periodic(v :: Vector, k :: Int)
    n = length(v)
    n % k != 0 && return false
    return all(i -> v[i] == v[mod(i+k,1:n)], 1 : n)
end

@doc raw"""
    period(v :: Vector)

Return the smallest positive integer `k` such that `v` is `k`-periodic.

"""
function period(v :: Vector)
    n = length(v)
    res = Int[]
    for k = 1 : n
        is_periodic(v, k) && return k
    end
end


@doc raw"""
    ehrhart_quasipolynomial_with_periods(P :: RationalPolygon{T}) where {T <: Integer}

Return a k x 3-matrix of normalized coefficients of the Ehrhart
quasipolynomial of a `k`-rational polygon `P`, together with a vector of it's
three periods.

"""
function ehrhart_quasipolynomial_with_periods(P :: RationalPolygon{T}) where {T <: Integer}
    k = rationality(P)
    M = Matrix{T}(undef, k, 3)
    A = normalized_area(P)
    ehrhart_values = [number_of_k_rational_points(P,t) for t = 1 : 2k]
    
    for t = 1 : k
        M[t,1] = A
        M[t,2] = -(2t+k)*A + 2k*(ehrhart_values[t+k] - ehrhart_values[t])
        M[t,3] = (t^2+t*k)*A + 2k*((t+k)*ehrhart_values[t] - t*ehrhart_values[t+k])
    end
    periods = [period(M[:,i]) for i = 1 : 3]

    return (M, periods)
end


@doc raw"""
    ehrhart_quasipolynomial(P :: RationalPolygon)

Return a k x 3-matrix of normalized coefficients of the Ehrhart quasipolynomial of a
`k`-rational polygon `P`.

"""
ehrhart_quasipolynomial(P :: RationalPolygon) =
ehrhart_quasipolynomial_with_periods(P)[1]


@doc raw"""
    ehrhart_quasipolynomial_periods(P :: RationalPolygon)

Return the periods of the Ehrhart quasipolynomial of a `k`-rational polygon
`P`.

"""
ehrhart_quasipolynomial_periods(P :: RationalPolygon) =
ehrhart_quasipolynomial_with_periods(P)[2]


@doc raw"""
    ehrhart_quasipolynomial_period(P :: RationalPolygon)

Return the period of the Ehrhart quasipolynomial of a `k`-rational polygon
`P`.

"""
ehrhart_quasipolynomial_period(P :: RationalPolygon) =
lcm(ehrhart_quasipolynomial_periods(P))


@doc raw"""
    is_quasiintegral(P :: RationalPolygon)

Check whether a rational polygon is quasiintegral, i.e. has an Ehrhart
polynomial.

"""
is_quasiintegral(P :: RationalPolygon) =
ehrhart_quasipolynomial_period(P) == 1
