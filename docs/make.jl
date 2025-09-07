if length(ARGS) > 1
    println("Usage: julia make.jl <format>")
    println("Format can be 'html', 'pdf', or 'tex'.")
    exit(1)
end

format = isempty(ARGS) ? "html" : ARGS[1]

if format ∉ ["html", "pdf", "tex"]
    println("Usage: julia make.jl <format>")
    println("Format can be 'html', 'pdf', or 'tex'.")
    exit(1)
end

@info "make.jl: Building documentation for format \"$format\""

using Documenter, DocumenterCitations, RationalPolygons, StaticArrays

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

if format == "html"
    makedocs(
        sitename = "RationalPolygons",
        pages = [
            "RationalPolygons.jl" => "index.md",
            "2D Geometry" => "2dgeometry.md",
            "Rational Polygons" => "polygons.md",
            "LDP polygons and toric surfaces" => "ldp.md",
            "Subpolygons" => "subpolygons.md",
            "Classifications" => "classifications.md",
            "Index" => "docs_index.md"
        ],
        plugins = [bib]
    )

    deploydocs(
        repo = "github.com/justus-springer/RationalPolygons.jl.git",
    )

elseif format == "pdf"
    makedocs(
        sitename = "RationalPolygons",
        format = Documenter.LaTeX(),
        pages = [
            "Home" => "index.md",
            "2D Geometry" => "2dgeometry.md",
            "Rational Polygons" => "polygons.md",
            "LDP polygons and toric surfaces" => "ldp.md",
            "Subpolygons" => "subpolygons.md",
            "Classifications" => "classifications.md",
            "Index" => "docs_index.md"
        ],
        plugins = [bib],
    )

elseif format == "tex"

else
    println("Invalid format. Please choose 'html', 'pdf', or 'tex'.")
    exit(1)
end
