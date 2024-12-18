@recipe function plot_recipe(Ps :: Union{RationalPolygon{T}, Vector{<:RationalPolygon{T}}}) where {T <: Integer}

    if Ps isa RationalPolygon
        Ps = [Ps]
    end

    k = rationality(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("Polygons must have the same rationality")

    xmin = minimum([v[1] for P ∈ Ps for v ∈ P]) - 1//k
    xmax = maximum([v[1] for P ∈ Ps for v ∈ P]) + 1//k
    ymin = minimum([v[2] for P ∈ Ps for v ∈ P]) - 1//k
    ymax = maximum([v[2] for P ∈ Ps for v ∈ P]) + 1//k

    lattice_points = [(x,y) for x = ceil(xmin) : floor(xmax) for y = ceil(ymin) : floor(ymax)]

    framestyle --> :none
    aspect_ratio --> true

    @series begin
        seriestype := :scatter
        markercolor --> :black
        xlims --> (xmin,xmax)
        ylims --> (ymin,ymax)
        xticks --> ceil(xmin):1:floor(xmax)
        yticks --> ceil(ymin):1:floor(ymax)
        aspect_ratio --> :equal
        axis --> false
        label --> false
        grid --> false
        lattice_points
    end

    for P ∈ Ps
        @series begin
            seriestype := :shape
            label --> false
            fillcolor --> :gray
            opacity --> 0.3
            [(v[1],v[2]) for v ∈ vertices(P)]
        end
    end

end


