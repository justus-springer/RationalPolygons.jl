@recipe function plot_recipe(Ps :: Union{RationalPolygon{T}, Vector{<:RationalPolygon{T}}}) where {T <: Integer}

    single = Ps isa RationalPolygon
    if single
        Ps = [Ps]
    end

    show_vertices = pop!(plotattributes, :show_vertices, false)

    k = lcm(rationality.(Ps))

    xmin = floor(minimum(v[1] for P ∈ Ps for v ∈ P)) - 1
    xmax = ceil(maximum(v[1] for P ∈ Ps for v ∈ P)) + 1
    ymin = floor(minimum(v[2] for P ∈ Ps for v ∈ P)) - 1
    ymax = ceil(maximum(v[2] for P ∈ Ps for v ∈ P)) + 1

    framestyle --> :none
    aspect_ratio --> :equal
    legend --> false

    for (i, P) ∈ enumerate(Ps)
        vs = [(v[1], v[2]) for v ∈ vertices(P)]
        color = single ? :gray : i
        @series begin
            seriestype := :shape
            fillcolor --> color
            fillalpha --> 0.3
            linecolor --> :black
            linewidth --> 1.5
            vs
        end
        if show_vertices
            @series begin
                seriestype := :scatter
                markercolor --> :white
                markerstrokecolor --> (single ? :black : color)
                markersize --> 6
                vs
            end
        end
    end

    if k > 1
        @series begin
            seriestype := :scatter
            markercolor --> :gray
            markerstrokewidth --> 0
            markersize --> 1.5
            [(x//k, y//k) for x = k*xmin : k*xmax for y = k*ymin : k*ymax if x % k != 0 || y % k != 0]
        end
    end

    @series begin
        seriestype := :scatter
        markercolor --> :black
        markersize --> 3
        [(x, y) for x = xmin : xmax for y = ymin : ymax]
    end

end
