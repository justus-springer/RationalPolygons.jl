module RationalPolygons

include("imports.jl")
include("exports.jl")

include("Point.jl")
include("graham.jl")
include("hilbert_basis.jl")

include("Line/Line.jl")
include("Line/intersect.jl")
include("Line/points_on_line_segment.jl")

include("AffineHalfplane/AffineHalfplane.jl")
include("AffineHalfplane/intersect.jl")

include("RationalPolygon/RationalPolygon.jl")
include("RationalPolygon/properties.jl")
include("RationalPolygon/io.jl")
include("RationalPolygon/h5.jl")
include("RationalPolygon/slice.jl")
include("RationalPolygon/lattice_width_data.jl")
include("RationalPolygon/interior_points.jl")
include("RationalPolygon/ehrhart.jl")
include("RationalPolygon/PolygonAutomorphismGroup.jl")
include("RationalPolygon/normal_form.jl")

include("Subpolygons/remove_vertex.jl")
include("Subpolygons/SubpolygonStorage.jl")
include("Subpolygons/InMemorySubpolygonStorage.jl")
include("Subpolygons/TextFilesSubpolygonStorage.jl")
include("Subpolygons/HDFSubpolygonStorage.jl")

include("Classification/one_interior_lattice_point.jl")
include("Classification/no_interior_lattice_points.jl")
include("Classification/lattice_polygons.jl")

end # module RationalPolygons
