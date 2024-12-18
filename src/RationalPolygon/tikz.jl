function tikz(io :: IO, Ps :: Vector{<:RationalPolygon{T}};
    gridline_width :: Tuple{Float64, Float64} = (1.0, 0.2),
    grid_size :: Tuple{Float64, Float64} = (0.5,0.5),
    line_width :: Float64 = 2.0,
    vertex_size :: Float64 = 2.0) where {T <: Integer}

    k = rationality(first(Ps))
    all(P -> rationality(P) == k, Ps) || error("Polygons must have the same rationality")

    xmin = minimum([k*v[1] for P ∈ Ps for v ∈ P]) - 0.5
    xmax = maximum([k*v[1] for P ∈ Ps for v ∈ P]) + 0.5
    ymin = minimum([k*v[2] for P ∈ Ps for v ∈ P]) - 0.5
    ymax = maximum([k*v[2] for P ∈ Ps for v ∈ P]) + 0.5

    println(io, "\\begin{tikzpicture}[x=$(grid_size[1])cm,y=$(grid_size[2])cm]\n")

    println(io, "% Gridlines ")
    println(io, "\\draw[step=$k,black,line width=$(gridline_width[1])pt] ($xmin,$ymin) grid ($xmax,$ymax);")
    println(io, "\\draw[help lines,step=1,black,line width=$(gridline_width[2])pt] ($xmin,$ymin) grid ($xmax,$ymax);")
    print(io, "\n\n")

    for i = 1 : length(Ps)

        P = Ps[i]
        N = number_of_vertices(P)
        V = vertex_matrix(P)

        println(io, "% Polygon no. $i")
        println(io, "\\fill[opacity=0.2]")
        for i = 1 : N
            print(io, "($(V[1,i]),$(V[2,i]))--")
        end
        println(io, "cycle;")

        println(io, "\\draw[line width=$(line_width)pt, color=black]")
        for i = 1 : N
            print(io, "($(V[1,i]),$(V[2,i]))--")
        end
        println(io, "cycle;")

        for i = 1 : N
            println(io, "\\draw [fill=black] ($(V[1,i]),$(V[2,i])) circle ($(vertex_size)pt);")
        end
        print(io, "\n\n")

    end

    println(io, "\\end{tikzpicture}")

end

tikz(io :: IO, P :: RationalPolygon{T}; kwargs...) where {T <: Integer} =
tikz(io, [P]; kwargs...)

