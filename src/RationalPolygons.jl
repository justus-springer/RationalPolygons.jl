module RationalPolygons

include("imports.jl")
include("exports.jl")

include("Point.jl")
include("Tools.jl")
include("hilbert_basis.jl")
include("Line.jl")
include("AffineHalfplane/AffineHalfplane.jl")
include("graham.jl")
include("RationalPolygon.jl")
include("io.jl")
include("slice.jl")
include("lattice_width_data.jl")
include("AffineHalfplane/intersect.jl")
include("interior_points.jl")
include("ehrhart.jl")
include("normal_form.jl")
include("PolygonAutomorphismGroup.jl")
include("Subpolygons/SubpolygonStorage.jl")
include("Subpolygons/subpolygons.jl")

include("Classification/classification.jl")
include("Classification/moving_out_the_edges.jl")

end # module RationalPolygons
