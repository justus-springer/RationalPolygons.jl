
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
        use_affine_normal_form :: Bool = false,
        hdf_path :: Union{Missing, String} = missing, 
        hdf_group :: String = "/",
        logging :: Bool = false) where {T <: Integer}

    if ismissing(hdf_path)
        st = InMemorySubpolygonStorage{T}(Ps; primitive, use_affine_normal_form)
        return subpolygons(st; logging)
    else
        f = h5open(hdf_path, "cw")
        if !haskey(f, hdf_group)
            create_group(f, hdf_group)
        end
        g = f[hdf_group]
        initialize_hdf_subpolygon_storage(g, Ps; primitive, use_affine_normal_form)
        return subpolygons(g; logging)
    end

end

subpolygons(P :: RationalPolygon{T},
        primitive :: Bool = false, 
        use_affine_normal_form :: Bool = false,
        hdf_path :: Union{Missing, String} = missing, 
        hdf_group :: String = "/",
        logging :: Bool = false) where {T <: Integer} =
subpolygons([P]; primitive, use_affine_normal_form, hdf_path, hdf_group, logging)
