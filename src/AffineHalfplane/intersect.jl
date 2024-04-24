
function _intersect_sorted_halfplanes(halfplanes :: Vector{<:AffineHalfplane{T}}) where {T <: Integer}

    result = empty(halfplanes)

    for H ∈ halfplanes
        while length(result) > 1 && intersection_point(line(result[end]), line(H)) ∉ result[end-1]
            pop!(result)
        end
        push!(result, H)
    end

    return result

end

function intersect_halfplanes(halfplanes :: Vector{<:AffineHalfplane{T}}) where {T <: Integer}

    n = length(halfplanes)
    sort!(halfplanes)

    upper_halfplanes, lower_halfplanes = AffineHalfplane{T}[], AffineHalfplane{T}[]
    for i = 1 : n
        H1, H2 = halfplanes[i], halfplanes[mod(i+1, 1:n)]
        a1, a2 = pseudo_angle(H1), pseudo_angle(H2)
        if a1 != a2
            push!(a1 ≥ 0 ? lower_halfplanes : upper_halfplanes, H1)
        end
    end

    upper_halfplanes = _intersect_sorted_halfplanes(upper_halfplanes)
    lower_halfplanes = _intersect_sorted_halfplanes(lower_halfplanes)

    found_right = false
    for i = 1 : length(upper_halfplanes), j = length(lower_halfplanes) : -1 : 1
        Hu, Hl = upper_halfplanes[i], lower_halfplanes[j]
        int = intersection_behaviour(line(Hu), line(Hl))
        int isa IntersectInPoint || continue
        if all(H -> int.p ∈ H, upper_halfplanes) && all(H -> int.p ∈ H, lower_halfplanes)
            found_right = true
            upper_halfplanes = upper_halfplanes[i : end]
            lower_halfplanes = lower_halfplanes[1 : j]
            break
        end
    end
    found_right || return EmptyPolygon{T}()

    found_left = false
    for i = length(upper_halfplanes) : -1 : 1, j = 1 : length(lower_halfplanes)
        Hu, Hl = upper_halfplanes[i], lower_halfplanes[j]
        int = intersection_behaviour(line(Hu), line(Hl))
        int isa IntersectInPoint || continue
        if all(H -> int.p ∈ H, upper_halfplanes) && all(H -> int.p ∈ H, lower_halfplanes)
            found_left = true
            upper_halfplanes = upper_halfplanes[1 : i]
            lower_halfplanes = lower_halfplanes[j : end]
            break
        end
    end
    found_left || return EmptyPolygon{T}()

    return IntersectionOfHalfplanes([upper_halfplanes ; lower_halfplanes])

end


