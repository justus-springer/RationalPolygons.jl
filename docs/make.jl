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

DocMeta.setdocmeta!(RationalPolygons, :DocTestSetup, :(using RationalPolygons, StaticArrays); recursive=true)

if format == "html"
    makedocs(
        modules = [RationalPolygons],
        sitename = "RationalPolygons",
        pages = [
            "RationalPolygons.jl" => "index.md",
            "2D Geometry" => "2dgeometry.md",
            "Rational Polygons" => "polygons.md",
            "Subpolygons" => "subpolygons.md",
            "Classifications" => "classifications.md",
            "Bibliography" => "bibliography.md",
            "Index" => "docs_index.md"
        ],
        warnonly = [:missing_docs, :doctest],
        plugins = [bib]
    )

    deploydocs(
        repo = "github.com/justus-springer/RationalPolygons.jl.git",
    )

elseif format == "thesis"
    makedocs(
        # modules = [RationalPolygons],
        sitename = "RationalPolygons",
        format = Documenter.LaTeX(platform = "none"),
        pages = [
            "Home" => "index.md",
            "2D Geometry" => "2dgeometry.md",
            "Rational Polygons" => "polygons.md",
            "Subpolygons" => "subpolygons.md",
            "Classifications" => "classifications.md",
        ],
        warnonly = :missing_docs,
        plugins = [bib],
    )

    filename = joinpath(@__DIR__, "build", "RationalPolygons.tex")

    txt = read(filename, String)

    # Add chapter heading
    txt = "\\chapter{RationalPolygons.jl}\n\\label{apx:julia_rational_polygons}\n" * txt

    # Add Tex root directive
    txt = "%!TEX root = thesis.tex\n\n" * txt

    # Fix displaying of emojis
    txt = replace(txt, "🗂️" => "|\\folder|", "🔢" => "|\\dataset|")

    # Fix citations
    txt = replace(txt, r"\[\\hyperref\[doc:(\w+)\]\{\d+\}\]" => s"\\cite{\1}")
    
    # Remvoe line breaks before equations
    txt = replace(txt, r"\n+(\\begin\{equation\*\})" => s"\n\1")

    # Remove math mode for references and add tilde
    txt = replace(txt, r" \\\( (\\ref\{\S+\}) \\\)" => s"~\1")

    write(filename, txt)

else
    println("Invalid format. Please choose 'html', 'pdf', or 'tex'.")
    exit(1)
end
