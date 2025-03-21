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
        "GKM varieties" => ["GKM" => "man/GKM/GKM.md", 
                            "Constructors" => "man/GKM/Constructors.md", 
                            "Properties" => "man/GKM/Properties.md", 
                            "Standard Constructions" => "man/GKM/STDconstructions.md",
                            "Operators" => "man/GKM/Operators.md", 
                            "Connections" => "man/GKM/Connections.md"],
        "GW invariants" => "man/GW.md",
        "References" => "references.md"]

makedocs(
    sitename = "GKMtest",
    format = Documenter.HTML(),
    modules = [GKMtest],
    warnonly = true,
    pages = pages,
    plugins=[bib]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
