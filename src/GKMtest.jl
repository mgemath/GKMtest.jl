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
include("flag.jl")
include("toric.jl")
include("betti.jl")
include("GKMsubgraphs.jl")
include("product.jl")
include("blowup.jl")
include("curveClasses.jl")

include("GW/includes.jl")
include("curveClasses.jl")

# files in the definitive format
# include("standard_constructions.jl")
# include("betti_numbers.jl")

end # module GKMtest
