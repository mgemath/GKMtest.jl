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
  con::Dict{Edge, Dict{Edge, Edge}} # assigns to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  a::Dict{Edge,Dict{Edge,ZZRingElem}} # w[e'_i] = w [e_i] - a_i * w[e]

  function GKM_connection(
    gkm::AbstractGKM_graph,
    con::Dict{Edge, Dict{Edge, Edge}},
    a::Dict{Edge,Dict{Edge,ZZRingElem}}
  )
    return new(gkm, con, a)
  end
end