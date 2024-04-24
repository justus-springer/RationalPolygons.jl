
struct EmptyPolygon{T <: Integer} <: RationalPolygon{T} end

Base.show(io :: IO, P :: EmptyPolygon) =
Base.print(io, "Empty polygon")

vertices(P :: EmptyPolygon{T}) where {T <: Integer} = RationalPoint{T}[]

edges(P :: EmptyPolygon{T}) where {T <: Integer} = Tuple{RationalPoint{T}, RationalPoint{T}}[]

affine_halfplanes(P :: EmptyPolygon{T}) where {T <: Integer} =
[AffineHalfplaneByNormalVector((one(T),zero(T)), one(T))
 AffineHalfplaneByNormalVector((-one(T),zero(T)), one(T))]

Base.in(x :: Point{T}, P :: EmptyPolygon) where {T <: Integer} = false
 
