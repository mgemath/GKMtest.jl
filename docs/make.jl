using Documenter, DocumenterCitations
using Oscar, GKMtest


bib = CitationBibliography(
    joinpath(@__DIR__, "src", "gkm_references.bib");
    style=:alpha
)

DocMeta.setdocmeta!(GKMtest, :DocTestSetup, :(using Oscar, GKMtest); recursive=true)


pages = [
        "Home" => "index.md",
        "GKM varieties" => ["GKM Graphs" => "GKM/GKM.md", 
                            "Constructors" => "GKM/Constructors.md", 
                            "Properties" => "GKM/Properties.md",
                            "Connections" => "GKM/Connections.md", 
                            "Standard Constructions" => "GKM/STDconstructions.md",
                            "Operators" => "GKM/Operators.md", 
                            "Cohomology" => "GKM/Cohomology.md",
                            "Curve Classes" => "GKM/CurveClasses.md",
                            "Vector Bundles" => "GKM/Vectorbundles.md",
                            "Seidel Space" => "GKM/Seidelspace.md"],
        "Gromov--Witten theory & Quantum Cohomology" => ["Gromov--Witten invariants" => "GW/GW.md",
                                                "Quantum Cohomology" => "GW/QH.md",
                                                "Seidel Elements / Shift Operators" => "GW/SeidelElements.md"],
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
    repo = "github.com/mgemath/GKMtest.jl.git",
    devbranch = nothing,
    devurl = "dev",
    versions = versions = ["stable" => "v^", "v#.#"],
)
