module RationalPolygons

include("imports.jl")
include("exports.jl")

include("Point.jl")
include("Tools.jl")
include("hilbert_basis.jl")
include("Line.jl")
include("AffineHalfplane/AffineHalfplane.jl")
include("graham.jl")
include("RationalPolygon/RationalPolygon.jl")
include("RationalPolygon/EmptyPolygon.jl")
include("RationalPolygon/ConvexHull.jl")
include("RationalPolygon/IntersectionOfHalfplanes.jl")
include("AffineHalfplane/intersect.jl")
include("interior_points.jl")
include("normal_form.jl")
include("subpolygons.jl")

include("Classification/one_interior_lattice_point.jl")

end # module RationalPolygons
