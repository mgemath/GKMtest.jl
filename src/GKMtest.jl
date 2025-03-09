@doc raw"""
**GKMtest** is a test module for GKM varieties to be included in OSCAR.

"""

module GKMtest

using Oscar, Combinatorics, ProgressMeter

include("exports.jl")
include("Types.jl")
include("GKMgraphs.jl")
include("GKMconnections.jl")
include("cohomology.jl")
# include("flag.jl") #obsolate: erase!
# include("toric.jl") #obsolate: erase!
include("betti.jl")
include("GKMsubgraphs.jl")
include("product.jl")
include("blowup.jl")
include("curveClasses.jl")

include("GP.jl")
include("different_w_types.jl")
include("GW/includes.jl")

# files in the definitive format
include("standard_constructions.jl")
# include("betti_numbers.jl")

end # module GKMtest
