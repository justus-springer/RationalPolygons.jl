
function subpolygons(
        st :: SubpolygonStorage{T};
        logging :: Bool = false) where {T <: Integer}

    while !is_finished(st)
        subpolygons_single_step(st; logging)
    end

    logging && @info "Found a total of $(total_count(st)) subpolygons"

    return return_value(st)

end


@doc raw"""
    subpolygons(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    subpolygons(P :: RationalPolygon{T}) where {T <: Integer}

Given a list of rational polygons with shared rationality and number of
interior lattice points, compute all subpolygons with the same number of
)nterior lattice points. The following keyword arguments are supported:

- `primitive :: Bool`: If set to true, only subpolygons with primitive vertices
are returned.

- `normal_form :: Symbol`: Can be either `:unimodular` or `:affine`. Used to
control which normal form should be used when comparing subpolygons for
equivalence.

- `out_path :: Union{Missing, String}`: If out_path is `missing`, then all
polygons are kept in memory. By specifying `out_path` to be a path to an empty
directory, the storage of the resulting polygons is delegated to the disk. In
this case, `subpolygons` will save the polygons to text files according to
their normalized area. Additionally, a file `last_volume` will be created and
constantly updated to hold the last fully completed area. This is useful for
resuming a lengthy computation at a later point, see also
`OnDiskSubpolygonStorage` for details on how to resume an unfinished
computations of subpolygons.

- `logging :: Bool`: Controls whether to show log messages showing the
progress.

"""
function subpolygons(Ps :: Vector{<:RationalPolygon{T}}; 
        primitive :: Bool = false, 
        normal_form :: Symbol = :unimodular,
        out_path :: Union{Missing, String} = missing, 
        logging :: Bool = false) where {T <: Integer}

    if ismissing(out_path)
        st = InMemorySubpolygonStorage{T}(Ps; primitive, normal_form)
    else
        st = OnDiskSubpolygonStorage{T}(Ps, out_path; primitive, normal_form)
    end

    return subpolygons(st; logging)
                    
end

subpolygons(P :: RationalPolygon{T},
        primitive :: Bool = false, 
        normal_form :: Symbol = :unimodular,
        out_path :: Union{Missing, String} = missing, 
        logging :: Bool = false) where {T <: Integer} =
subpolygons([P]; primitive, normal_form, out_path, logging)
