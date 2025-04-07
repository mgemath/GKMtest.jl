@doc raw"""
**GKMtest** is a test module for GKM varieties to be included in OSCAR.

"""

module GKMtest

using Oscar, Combinatorics, ProgressMeter

## GKM
include("exports.jl")
include("imports.jl")
include("Types.jl")
include("different_w_types.jl")

## Constructors 
## Properties
include("GKMgraphs.jl")
include("betti.jl")

## Standard Constructions
include("standard_constructions.jl")
include("GP.jl")

## Operators
include("GKMsubgraphs.jl")
include("product.jl")
include("blowup.jl")

## Connections
include("GKMconnections.jl")

## Cohomology
include("cohomology.jl")
include("curveClasses.jl")

## Vector Bundles
include("equivariant_bundles.jl")

## Seidel Space
include("Seidel_space.jl")


## GW
include("GW/includes.jl")

## obsolate
include("obsolate/obsolate.jl")

## experimental
include("bruhat.jl")


end # module GKMtest
