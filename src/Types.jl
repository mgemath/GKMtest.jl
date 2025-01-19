###############################################################################
#
#   GKMtest
#
###############################################################################


@attributes mutable struct AbstractGKM_graph
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

@attributes mutable struct AbstractGKM_subgraph
  super::AbstractGKM_graph
  self::AbstractGKM_graph # the GKM subgraph which forgets about the supergraph
  vDict::Vector{Int64} # track how vertices of the subgraph are mapped to that of the supergraph (since Oscar always uses {1, ..., n} as vertex set)

  function AbstractGKM_subgraph(super::AbstractGKM_graph, self::AbstractGKM_graph, vDict::Vector{Int64})
    return new(super, self, vDict)
  end
end

struct GKM_cohomology_ring
  gkm::AbstractGKM_graph
  coeffRing::QQMPolyRing # H_T^*(point;Q)
  cohomRing::FreeMod{QQMPolyRingElem} # H_T^*(X; Q), but without checks for consistency (see isGKMclass in cohomology.jl)

  function GKM_cohomology_ring(
    gkm::AbstractGKM_graph,
    coeffRing::QQMPolyRing,
    cohomRing::FreeMod{QQMPolyRingElem}
  )
    return new(gkm, coeffRing, cohomRing)
  end
end

mutable struct GKM_connection
  gkm::AbstractGKM_graph
  con::Dict{Tuple{Edge, Edge}, Edge} # assigns to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} # w[e'_i] = w [e_i] - a_i * w[e]

  function GKM_connection(
    gkm::AbstractGKM_graph,
    con::Dict{Tuple{Edge, Edge}, Edge},
    a::Dict{Tuple{Edge, Edge},ZZRingElem}
  )
    return new(gkm, con, a)
  end
end