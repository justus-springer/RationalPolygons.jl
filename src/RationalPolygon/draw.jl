function plot_polygon(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}
    k = rationality(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("Polygons must have the same rationality")

    xmin = minimum([v[1] for P ∈ Ps for v ∈ P]) - 1//k
    xmax = maximum([v[1] for P ∈ Ps for v ∈ P]) + 1//k
    ymin = minimum([v[2] for P ∈ Ps for v ∈ P]) - 1//k
    ymax = maximum([v[2] for P ∈ Ps for v ∈ P]) + 1//k

    lattice_points = [(x,y) for x = ceil(xmin) : floor(xmax) for y = ceil(ymin) : floor(ymax)]

    plot(lattice_points,
         seriestype = :scatter,
         markercolor = :black,
         xlims = (xmin,xmax),
         ylims = (ymin,ymax),
         xticks = ceil(xmin):1:floor(xmax),
         yticks = ceil(ymin):1:floor(ymax),
         aspect_ratio = :equal,
         axis = false,
         label = false,
         grid = false)

    P = first(Ps)
    for P ∈ Ps
        vs = vertices(P)
        k = rationality(P)
        shape = Shape([(v[1],v[2]) for v ∈ vs])
        plot!(shape, label = false, fillcolor = plot_color(:gray, 0.3), show = true)
    end

end

plot_polygon(P :: RationalPolygon{T}) where {T <: Integer} = plot_polygon([P])
