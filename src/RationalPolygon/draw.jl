function plot_polygon(P :: RationalPolygon{T,N}) where {N, T <: Integer}
    vs = vertices(P)
    k = rationality(P)
    shape = Shape([(v[1],v[2]) for v ∈ vs])
    xmin = minimum([v[1] for v ∈ vs]) - 1//k
    xmax = maximum([v[1] for v ∈ vs]) + 1//k
    ymin = minimum([v[2] for v ∈ vs]) - 1//k
    ymax = maximum([v[2] for v ∈ vs]) + 1//k
    plot(shape, 
         xlims = (xmin,xmax),
         ylims = (ymin,ymax),
         xticks = ceil(xmin):1:floor(xmax),
         yticks = ceil(ymin):1:floor(ymax),
         aspect_ratio = :equal,
         axis = false,
         label = false,
         grid = false,
         fillcolor = plot_color(:gray, 0.3))

    lattice_points = [(x,y) for x = ceil(xmin) : floor(xmax) for y = ceil(ymin) : floor(ymax)]

    plot!(lattice_points,
          seriestype = :scatter,
          markercolor = :black,
          label = :false)
end

