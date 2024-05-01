
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

@doc raw"""
    is_bounded(halfplanes :: Vector{<:AffineHalfplane{T}}) where {T <: Integer}   

Check whether the intersection of the given halfplanes is bounded, i.e.
describes a polygon.

"""
function is_bounded(halfplanes :: Vector{<:AffineHalfplane{T}}) where {T <: Integer}
    angles = sort!(pseudo_angle.(halfplanes))
    for i = 1 : length(angles)-1
        a1, a2 = angles[i], angles[i+1]
        a2 - a1 >= 2 && return false
    end
    a1, a2 = last(angles), first(angles)
    (a2+4) - a1 >= 2 && return false
    return true
end

function intersect_halfplanes(halfplanes :: Vector{<:AffineHalfplane{T}}; rationality :: Union{Missing,T} = missing) where {T <: Integer}

    n = length(halfplanes)

    is_bounded(halfplanes) || error("the intersection of these halfplanes is unbounded")

    # sort the halfplanes, see the `isless` implementation of
    # `AffineHalfplanes`: Halfplanes are sorted first by their angle and
    # second by their distance from the origin
    sort!(halfplanes)

    # split up into lower and upper halfplanes. At the same time,
    # discard any halfplanes that superset of other halfplanes and thus,
    # redundant. Since halfplanes of the same angle are sorted by
    # distance, we only keep the halfplane with the largest distance
    # for every given angle.
    upper_halfplanes, lower_halfplanes = AffineHalfplane{T}[], AffineHalfplane{T}[]
    for i = 1 : n
        H1, H2 = halfplanes[i], halfplanes[mod(i+1, 1:n)]
        a1, a2 = pseudo_angle(H1), pseudo_angle(H2)
        if a1 ≠ a2
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

    new_halfplanes = [upper_halfplanes ; lower_halfplanes]
    r = length(new_halfplanes)

    # construct the vertices of the resulting polygon
    vs = RationalPoint{T}[]
    for i = 1 : r
        H1, H2 = new_halfplanes[i], new_halfplanes[mod(i+1,1:r)]
        push!(vs, intersection_point(line(H1),line(H2)))
    end
    if ismissing(rationality)
        P = convex_hull(vs)
    else
        P = convex_hull(vs, rationality)
    end
    set_attribute!(P, :affine_halfplanes, new_halfplanes)

    return P
    
end


