using Documenter, DocumenterCitations, RationalPolygons

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    sitename = "RationalPolygons",
    pages = [
        "Home" => "index.md",
        "2D Geometry" => "2dgeometry.md",
        "Rational Polygons" => "polygons.md",
        "Subpolygons" => "subpolygons.md",
        "Classifications" => "classifications.md",
        "Index" => "docs_index.md"
    ],
    plugins = [bib]
)

deploydocs(
    repo = "github.com/justus-springer/RationalPolygons.jl.git",
)