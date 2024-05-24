
export classify_maximal_lattice_polygons_with_collinear_interior_points
function classify_maximal_lattice_polygons_with_collinear_interior_points(g :: Int, T :: Type{<:Integer} = Int)
    Ps = RationalPolygon{T}[]
    push!(Ps,RationalPolygon(LatticePoint{T}[(0,-1),(2g+2,-1),(0,1)], 1))
    for i = 1 : g+1
        push!(Ps,RationalPolygon(LatticePoint{T}[(0,-1),(2g+2-i,-1),(i,1),(0,1)], 1))
    end
    return Ps
end

export classify_empty_lattice_polygons_by_boundary_points
function classify_empty_lattice_polygons_by_boundary_points(R :: Int, T :: Type{<:Integer} = Int)
    Ps = RationalPolygon{T}[]
    push!(Ps, RationalPolygon(LatticePoint{T}[(0,0),(R-2,0),(0,1)], 1))
    for i = 1 : (R÷2)-1
        push!(Ps, RationalPolygon(LatticePoint{T}[(0,0),(R-i-2,0),(i,1),(0,1)], 1))
    end
    R == 6 && push!(Ps, RationalPolygon(LatticePoint{T}[(0,0),(2,0),(0,2)], 1))
    return Ps
end

export classify_maximal_lattice_polygons_with_empty_fine_interior
function classify_maximal_lattice_polygons_with_empty_fine_interior(g :: Int, T :: Type{<:Integer} = Int)
    g == 2 && return RationalPolygon{T}[]
    Ps = classify_empty_lattice_polygons_by_boundary_points(g, T)
    Qs = move_out_edges.(Ps)
    filter!(Q -> rationality(Q) == 1, Qs)
    return Qs
end

export _add_to_polygons_dict
function _add_to_polygons_dict(d :: Dict{Int,Vector{RationalPolygon{T}}}, P :: RationalPolygon{T}) where {T <: Integer}
    i = number_of_interior_lattice_points(P)
    if !haskey(d, i)
        d[i] = RationalPolygon{T}[]
    end
    push!(d[i], P)
end

export classify_lattice_polygons_up_to_genus
function classify_lattice_polygons_up_to_genus(g :: Int, T :: Type{<:Integer} = Int)
    res = Dict{Int,Vector{RationalPolygon{T}}}()

    res[1] = classify_polygons_genus_one(1)
    for P ∈ res[1]
        Q = move_out_edges(P)
        _add_to_polygons_dict(res, Q)
    end

    for i = 2 : g
        @info i
        if !haskey(res,i)
            res[i] = RationalPolygon{T}[]
        end
        append!(res[i], classify_maximal_lattice_polygons_with_empty_fine_interior(i,T))
        append!(res[i], classify_maximal_lattice_polygons_with_collinear_interior_points(i,T))
        Ps = subpolygons(res[i]; normal_form = :affine)
        for P ∈ Ps
            Q = move_out_edges(P)
            if rationality(Q) == 1
                _add_to_polygons_dict(res, Q)
            end
        end
        res[i] = Ps
    end

    return res

end

