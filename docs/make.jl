if length(ARGS) > 1
    println("Usage: julia make.jl <format>")
    println("Format can be 'html', 'thesis'")
    exit(1)
end

format = isempty(ARGS) ? "html" : ARGS[1]

if format ∉ ["html", "thesis"]
    println("Usage: julia make.jl <format>")
    println("Format can be 'html', 'thesis'")
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
            "Bibliography" => "bibliography.md",
            "Index" => "docs_index.md"
        ],
        plugins = [bib]
    )

    deploydocs(
        repo = "github.com/justus-springer/RationalPolygons.jl.git",
    )

elseif format == "thesis"
    makedocs(
        sitename = "RationalPolygons",
        format = Documenter.LaTeX(platform = "none"),
        pages = [
            "2D Geometry" => "2dgeometry.md",
            "Rational Polygons" => "polygons.md",
            "LDP polygons and toric surfaces" => "ldp.md",
            "Subpolygons" => "subpolygons.md",
            "Classifications" => "classifications.md",
        ],
        plugins = [bib],
    )

    filename = joinpath(@__DIR__, "build", "RationalPolygons.tex")

    txt = read(filename, String)

    # Add preamble
    txt = "%!TEX root = thesis.tex\n" * txt

    # Fix displaying of emojis
    txt = replace(txt, "🗂️" => "|\\folder|", "🔢" => "|\\dataset|")

    # Fix citations
    txt = replace(txt, r"\[\\hyperref\[doc:(\w+)\]\{\d+\}\]" => s"\\cite{\1}")

    write(filename, txt)

else
    println("Invalid format. Please choose 'html', 'pdf', or 'tex'.")
    exit(1)
end
