
@doc raw"""
    k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}
    
Return all `k`-rational points lying on the line segment from p to q.
If `interior` is set to false, this includes `p` and `q` themselves, if
they are `k`-rational. Otherwise, `p` and `q` are always excluded.

"""
function k_rational_points_on_line_segment(k :: T, p :: Point{T}, q :: Point{T}; interior = true) where {T <: Integer}

    if q[1] < p[1] || (p[1] == q[1] && q[2] < p[2])
        return k_rational_points_on_line_segment(k, q, p; interior)
    end

    L = line_through_points(p,q)
    res = RationalPoint{T}[]
    if p[1] < q[1]
        lower_bound = interior ? floor(T, k*p[1]+1) // k : ceil(T, k*p[1]) // k
        upper_bound = interior ? ceil(T, k*q[1]-1) // k : floor(T, k*q[1]) // k
        for t = lower_bound : 1 // k : upper_bound
            x = intersection_point(L, vertical_line(t))
            is_k_rational(k,x) || continue
            push!(res,x)
        end
    elseif p[2] < q[2]
        lower_bound = interior ? floor(T, k*p[2]+1) // k : ceil(T, k*p[2]) // k
        upper_bound = interior ? ceil(T, k*q[2]-1) // k : floor(T, k*q[2]) // k
        for t = lower_bound : 1 // k : upper_bound
            x = intersection_point(L, horizontal_line(t))
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


