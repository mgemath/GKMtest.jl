"""
**GKMtest** is a test module for GKM varieties to be included in OSCAR.

"""

module GKMtest

using Oscar, Combinatorics

include("exports.jl")
include("Types.jl")
include("GKMgraphs.jl")
include("GKMconnections.jl")
include("cohomology.jl")

include("GW/includes.jl")

greet() = print("GKMtest loaded")

end # module GKMtest
