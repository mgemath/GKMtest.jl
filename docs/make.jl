using Documenter
using GKMtest
using DocumenterCitations

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
                            "Standard Constructions" => "GKM/STDconstructions.md",
                            "Operators" => "GKM/Operators.md", 
                            "Connections" => "GKM/Connections.md",
                            "Cohomology" => "GKM/Cohomology.md"],
        "GW invariants" => "GW/GW.md",
        "References" => "references.md"]

makedocs(
    sitename = "GKMtest",
    format = Documenter.HTML(),
    modules = [GKMtest],
    warnonly = true,
    pages = pages,
    plugins=[bib],
    doctest = false,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mgemath/GKMtest.jl.git"
)
