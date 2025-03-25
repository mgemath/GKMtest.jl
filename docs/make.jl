using Documenter, DocumenterCitations
using GKMtest


bib = CitationBibliography(
    joinpath(@__DIR__, "src", "gkm_references.bib");
    style=:alpha
)

DocMeta.setdocmeta!(GKMtest, :DocTestSetup, :(using Oscar, GKMtest); recursive=true)


pages = [
        "Home" => "index.md",
        "GKM varieties" => ["GKM" => "GKM/GKM.md", 
                            "Constructors" => "GKM/Constructors.md", 
                            "Properties" => "GKM/Properties.md",
                            "Connections" => "GKM/Connections.md", 
                            "Standard Constructions" => "GKM/STDconstructions.md",
                            "Operators" => "GKM/Operators.md", 
                            "Cohomology" => "GKM/Cohomology.md",
                            "Vector Bundles" => "GKM/Vectorbundles.md",
                            "Seidel Space" => "GKM/Seidelspace.md"],
        "GW invariants" => "GW/GW.md",
        "References" => "references.md"]

makedocs(
    sitename = "GKMtest",
    format = Documenter.HTML(),
    modules = [GKMtest],
    warnonly = true,
    pages = pages,
    plugins=[bib],
    doctest = true,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mgemath/GKMtest.jl.git"
)
