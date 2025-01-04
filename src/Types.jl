###############################################################################
#
#   GKMtest
#
###############################################################################

# GKM_graph{C<:RingElem} = Tuple{Graph{Directed}, Dict{Edge, AbstractAlgebra.Generic.FreeModule{C}}}

@attributes mutable struct AbstractGKM_graph #<: GKM_graph{C}
  g::Graph
  labels::Vector{String}
  M::AbstractAlgebra.Generic.FreeModule{ZZRingElem} # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}} # weight of the T-action

  function AbstractGKM_graph(
    g::Graph,
    labels::Vector{String},
    M::AbstractAlgebra.Generic.FreeModule{ZZRingElem}, # character group
    w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}
  )
    return new(g, labels, M, w)
  end
end