module RationalPolygons

include("imports.jl")
include("exports.jl")

include("Point.jl")
include("graham.jl")
include("hilbert_basis.jl")
include("snf.jl")
include("hnf.jl")

include("Line/Line.jl")
include("Line/intersect.jl")
include("Line/points_on_line_segment.jl")

include("AffineHalfplane/AffineHalfplane.jl")
include("AffineHalfplane/intersect.jl")

include("RationalPolygon/RationalPolygon.jl")
include("RationalPolygon/properties.jl")
include("RationalPolygon/ldp.jl")
include("RationalPolygon/io.jl")
include("RationalPolygon/h5.jl")
include("RationalPolygon/slice.jl")
include("RationalPolygon/lattice_width_data.jl")
include("RationalPolygon/interior_points.jl")
include("RationalPolygon/ehrhart.jl")
include("RationalPolygon/PolygonAutomorphismGroup.jl")
include("RationalPolygon/normal_form.jl")
include("RationalPolygon/plot_recipe.jl")
include("RationalPolygon/tikz.jl")

include("Subpolygons/remove_vertex.jl")
include("Subpolygons/SubpolygonStorage.jl")
include("Subpolygons/InMemorySubpolygonStorage.jl")
include("Subpolygons/HDFSubpolygonStorage.jl")

include("Classification/collinear_interior_points.jl")
include("Classification/one_interior_lattice_point.jl")
include("Classification/no_interior_lattice_points.jl")
include("Classification/castryck.jl")
include("Classification/koelman.jl")
include("Classification/baeuerle.jl")
include("Classification/picard_index_triangles.jl")
include("Classification/kasprzyk_kreuzer_nill.jl")
include("Classification/triangles_integral_degree.jl")
include("Classification/quadrilaterals.jl")

end # module RationalPolygons
