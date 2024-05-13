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
include("AffineHalfplane/intersect.jl")
include("interior_points.jl")
include("normal_form.jl")
include("Subpolygons/subpolygons.jl")
include("Subpolygons/SubpolygonStorage.jl")

include("classification.jl")

end # module RationalPolygons
