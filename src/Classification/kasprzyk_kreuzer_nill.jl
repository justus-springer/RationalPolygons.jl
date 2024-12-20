export classify_lattice_polygons_by_gorenstein_index
export choose_next_vertex
export initial_special_facets


struct PartialLDP{T<:Integer,N,M}
    P :: RationalPolygon{T,N,M}
    initial_local_index :: T
    ymin :: T
end

function initial_special_facets(index :: T, initial_local_index :: T) where {T <: Integer}
    res = PartialLDP{T,2,4}[]
    for a = 0 : initial_local_index - 1
        gcd(a, initial_local_index) == 1 || continue
        for b = -a+1 : 2*index*(initial_local_index+1) - a
            gcd(b, initial_local_index) == 1 || continue
            newP = RationalPolygon(SMatrix{2,2,T}(b, initial_local_index, -a, initial_local_index), 1)
            newLDP = PartialLDP{T,2,4}(newP, initial_local_index, -initial_local_index * (initial_local_index + 1))

            push!(res, newLDP)
        end
    end
    return res
end

function choose_next_vertex(ldps :: Vector{<:PartialLDP{T,N}}, index :: T) where {N, T <: Integer}

    res = PartialLDP{T,N+1,2*(N+1)}[]
    for ldp ∈ ldps

        P, local_index, ymin = ldp.P, ldp.initial_local_index, ldp.ymin

        V = vertex_matrix(P)
        b, a = V[1,N-1], -V[1,N]
        v1 = scaled_vertex(P,1)

        H1 = affine_halfplane(v1, LatticePoint{T}(0,0))
        H2 = affine_halfplane(P,1)
        H3 = -affine_halfplane(P,N)

        Hbot = affine_halfplane(RationalPoint{T}(0,1), ymin)
        Htop = affine_halfplane(RationalPoint{T}(0,-1), 1 - local_index)

        # Halfplanes from Cor. 6.2
        H4 = affine_halfplane(RationalPoint{T}(local_index, a - index), - index * local_index)
        H5 = affine_halfplane(RationalPoint{T}(-local_index, b - index), -index * local_index)

        search_space = intersect_halfplanes([Hbot, Htop, H1, H2, H3, H4, H5])

        for v ∈ lattice_points(search_space)
            is_primitive(v) || continue

            contains_in_interior(v, H1) || continue
            contains_in_interior(v, H2) || continue
            contains_in_interior(v, H3) || continue

            # The new local index must divide the given index
            new_local_index = det(v1,v) ÷ gcd(v[2] - v1[2], v1[1] - v[1])
            index % new_local_index == 0 || continue

            new_polygon = RationalPolygon([v V], 1)
            new_ymin = ymin - min(v[2], 0) + sum([j for j = max(min(v[2], v1[2]),0) + 1 : max(v[2], v1[2]) -1])
            new_ldp = PartialLDP{T,N+1,2*(N+1)}(new_polygon, local_index, new_ymin)

            push!(res, new_ldp)
        end
    end

    return res

end

function classify_lattice_polygons_by_gorenstein_index(index :: T) where {T <: Integer}
    Pss = Set{<:RationalPolygon{T}}[]

    for local_index = 1 : index
        index % local_index == 0 || continue

        initial_facets = initial_special_facets(index, local_index)
        ldps = choose_next_vertex(initial_facets, index)
        N = 3
        while !isempty(ldps)
            @info local_index, N, length(ldps)
            length(Pss) < N-2 && push!(Pss, Set{RationalPolygon{T,N,2*N}}())

            # Save LDP polygons from previous interation
            new_ldps = [unimodular_normal_form(ldp.P) for ldp ∈ ldps if is_ldp(ldp.P) && gorenstein_index(ldp.P) == index]
            union!(Pss[N-2], new_ldps)

            # Choose next vertices
            ldps = choose_next_vertex(ldps, index)
            N += 1
        end

    end

    return Pss

end



