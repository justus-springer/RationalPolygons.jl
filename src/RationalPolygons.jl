module RationalPolygons

include("imports.jl")
include("exports.jl")

include("Point.jl")
include("Tools.jl")
include("Line.jl")
include("AffineHalfplane/AffineHalfplane.jl")
include("graham.jl")
include("RationalPolygon/RationalPolygon.jl")
include("RationalPolygon/ConvexHull.jl")
include("RationalPolygon/IntersectionOfHalfplanes.jl")
include("RationalPolygon/EmptyPolygon.jl")
include("AffineHalfplane/intersect.jl")
include("interior_points.jl")

end # module RationalPolygons
