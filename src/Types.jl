###############################################################################
#
#   GKMtest
#
###############################################################################

# We introduce abstract types so that we can 
abstract type AbstractGKM_H2 end
abstract type AbstractGKM_connection end
abstract type AbstractGKM_cohomology_ring end

@attributes mutable struct AbstractGKM_graph
  g::Graph
  labels::Vector{String}
  M::AbstractAlgebra.Generic.FreeModule{ZZRingElem} # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}} # weight of the T-action
  # This should always be set and can be accessed directly:
  # It should not be changed.
  equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring} # actual type will be Union{Nothing, GKM_cohomology_ring}
  # Use GKM_second_homology() to access this:
  curveClasses::Union{Nothing, AbstractGKM_H2} # actual type will be Union{Nothing, GKM_H2}
  # Use get_GKM_connection() to access this:
  connection::Union{Nothing, AbstractGKM_connection} # actual type will be Union{Nothing, GKM_connection}

  function AbstractGKM_graph(
    g::Graph,
    labels::Vector{String},
    M::AbstractAlgebra.Generic.FreeModule{ZZRingElem}, # character group
    w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}},
    equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring},
    curveClasses::Union{Nothing, AbstractGKM_H2},
    connection::Union{Nothing, AbstractGKM_connection}
  )
    return new(g, labels, M, w, equivariantCohomology, curveClasses, connection)
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

struct GKM_cohomology_ring <: AbstractGKM_cohomology_ring
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

mutable struct GKM_connection <: AbstractGKM_connection
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

mutable struct GKM_H2 <: AbstractGKM_H2
  gkm::AbstractGKM_graph
  edgeLattice::AbstractAlgebra.FPModule{ZZRingElem}
  H2::AbstractAlgebra.FPModule{ZZRingElem} # quotient of edgeLattice by relations in H_2.
  edgeToGenIndex::Dict{Edge, Int64}
  quotientMap::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem} # Z-module homomorphism from edgeLattice to H2
  dualConeRaySum::RayVector{QQFieldElem} # sum of rays of the dual cone of the edgeCurveClasses, normalized so that the minimum of pairings with edge curve classes is 1.
  dualCone::Cone{QQFieldElem} # dual cone of the cone of effective curve classes
  chernNumber::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem} # Z-module homomorphism from H2 to ZZ, giving the curve class evaluated on the first chern class of the tangent bundle of the space.

  function GKM_H2(
    gkm::AbstractGKM_graph,
    edgeLattice::AbstractAlgebra.FPModule{ZZRingElem},
    H2::AbstractAlgebra.FPModule{ZZRingElem},
    edgeToGenIndex::Dict{Edge, Int64},
    quotientMap::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem},
    dualConeRaySum::RayVector{QQFieldElem},
    dualCone::Cone{QQFieldElem},
    chernNumber::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}
  )
    return new(gkm, edgeLattice, H2, edgeToGenIndex, quotientMap, dualConeRaySum, dualCone, chernNumber)
  end
end

CurveClass_type = AbstractAlgebra.FPModuleElem{ZZRingElem}