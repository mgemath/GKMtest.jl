using Documenter
using GKMtest

DocMeta.setdocmeta!(GKMtest, :DocTestSetup, :(using Oscar, GKMtest); recursive=true)

makedocs(
    sitename = "GKMtest",
    format = Documenter.HTML(),
    modules = [GKMtest],
    warnonly = true,
    pages = [
        "Home" => "index.md",
        "GKM varieties" => "man/GKM.md",
        "GW invariants" => "man/GW.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
