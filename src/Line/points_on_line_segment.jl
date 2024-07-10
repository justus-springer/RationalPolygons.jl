
@doc raw"""
    k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}
    
Return all `k`-rational points lying on the line segment from p to q.
If `interior` is set to false, this includes `p` and `q` themselves, if
they are `k`-rational. Otherwise, `p` and `q` are always excluded.

"""
function k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}

    if p == q
        if interior || !is_k_rational(k,p)
            return RationalPoint{T}[]
        else
            return RationalPoint{T}[p]
        end
    end

    if q[1] < p[1] || (p[1] == q[1] && q[2] < p[2])
        return k_rational_points_on_line_segment(k, q, p; interior)
    end

    L = line_through_points(p,q)
    res = RationalPoint{T}[]
    if p[1] < q[1]
        lower_bound = interior ? floor(T, k*p[1]+1) // k : ceil(T, k*p[1]) // k
        upper_bound = interior ? ceil(T, k*q[1]-1) // k : floor(T, k*q[1]) // k
        for t = lower_bound : 1 // k : upper_bound
            x = intersection_point(L, vertical_line(Rational{T}(t)))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    elseif p[2] < q[2]
        lower_bound = interior ? floor(T, k*p[2]+1) // k : ceil(T, k*p[2]) // k
        upper_bound = interior ? ceil(T, k*q[2]-1) // k : floor(T, k*q[2]) // k
        for t = lower_bound : 1 // k : upper_bound
            x = intersection_point(L, horizontal_line(Rational{T}(t)))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    end

    return res

end


@doc raw"""
    integral_points_on_line_segment(p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}

Return all integral points lying on the line segment from p to q. If `interior`
is set to false, this includes `p` and `q` themselves, if they are integral.
Otherwise, `p` and `q` are always excluded.

"""
integral_points_on_line_segment(p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer} =
k_rational_points_on_line_segment(1, p, q; interior)


@doc raw"""
    integral_points_on_line_segment_with_given_integral_point(p :: Point{T}, q :: Point{T}, x0 :: LatticePoint{T})

Given a lattice point `x0` lying on the line segment between the points `p` and
`q`, return all lattice points on that line segment.

"""
function integral_points_on_line_segment_with_given_integral_point(p :: Point{T}, q :: Point{T}, x0 :: LatticePoint{T}) where {T <: Integer}

    if q[1] < p[1] || (p[1] == q[1] && q[2] < p[2])
        return integral_points_on_line_segment_with_given_integral_point(q, p,x0)
    end

    p == q && return [x0]

    v = primitivize(q - p)
    if v[1] != 0
        lower_bound = cld(p[1] - x0[1], v[1])
        upper_bound = fld(q[1] - x0[1], v[1])
        return [x0 + i*v for i = lower_bound : upper_bound]
    elseif v[2] != 0
        lower_bound = cld(p[2] - x0[2], v[2])
        upper_bound = fld(q[2] - x0[2], v[2])
        return [x0 + i*v for i = lower_bound : upper_bound]
    end

end
