
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

    # sort the halfplanes, see the `isless` implementation of
    # `AffineHalfplanes`: Halfplanes are sorted first by their angle and
    # second by their distance from the origin
    sort!(halfplanes)

    # boundedness check
    # There are other ways in which the given halfplanes do not give
    # a pounded polytope, for example if there are two parallel horizontal
    # halfplanes. This is checked for later, when intersecting the upper
    # with the lower halfplanes
    !isempty(filter(H -> -2 < pseudo_angle(H) < 0, halfplanes)) || error("this polyhedra is unbounded")
    !isempty(filter(H -> 0 < pseudo_angle(H) < 2, halfplanes)) || error("this polyhedra is unbounded")

    # split up into lower and upper halfplanes. At the same time,
    # discard any halfplanes that superset of other halfplanes and thus,
    # redundant. Since halfplanes of the same angle are sorted by
    # distance, we only keep the halfplane with the largest distance
    # for every given angle.
    upper_halfplanes, lower_halfplanes = AffineHalfplane{T}[], AffineHalfplane{T}[]
    for i = 1 : n
        H1, H2 = halfplanes[i], halfplanes[mod(i+1, 1:n)]
        a1, a2 = pseudo_angle(H1), pseudo_angle(H2)
        if a1 != a2
            push!(a1 ≥ 0 ? lower_halfplanes : upper_halfplanes, H1)
        end
    end

    # remove redundancies in upper and lower halfplanes respectively.
    # This is basically the dual the graham scan algorithm.
    upper_halfplanes = _intersect_sorted_halfplanes(upper_halfplanes)
    lower_halfplanes = _intersect_sorted_halfplanes(lower_halfplanes)


    # Now we need to intersect the upper halfplanes with the lower
    # halfplanes. For this, we determine intesection points of the
    # lines and check if it lies in all the halfplanes.
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


