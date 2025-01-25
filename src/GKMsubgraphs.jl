import Oscar.has_edge
import Oscar.has_vertex

"""
Return the GKM subgraph induced by the given vertices
This does not check if the result is a valid GKM graph.
"""
function GKMsubgraph_from_vertices(gkm::AbstractGKM_graph, vertices::Vector{Int64}) :: AbstractGKM_subgraph

  vertices = sort(vertices)
  vDict = zeros(Int64, 0)

  last = 0
  for v in vertices
    @req v > 0 "Vertex index must be positive"
    if v == last #don't count duplicate vertices
      continue
    end
    push!(vDict, v)
    last = v
  end

  subnv = length(vDict)
  labels = [gkm.labels[vDict[i]] for i in 1:subnv]
  subGKM = gkm_graph(Graph{Undirected}(subnv), labels, gkm.M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()) # use the same character lattice
  
  for e in edges(gkm.g)
    if src(e) in vDict && dst(e) in vDict
      sd = indexin([src(e), dst(e)], vDict)
      GKMadd_edge!(subGKM, sd[1], sd[2], gkm.w[e])
    end
  end
  return AbstractGKM_subgraph(gkm, subGKM, vDict)
end

"""
Return the GKM subgraph induced by the given vertex labels
This does not check if the result is a valid GKM graph.
"""
function GKMsubgraph_from_vertices(gkm::AbstractGKM_graph, vertexLabels::Vector{String}) :: AbstractGKM_subgraph
  for l in vertexLabels
    @req l in gkm.labels "Label $l not found"
  end
  vertices::Vector{Int64} = indexin(vertexLabels, gkm.labels) # need to specify Vector{Int64} as indexin returns vector of Union{Nothing, Int64}.
  return GKMsubgraph_from_vertices(gkm, vertices)
end

"""
Return the GKM subgraph induced by the given edges
This does not check if the result is a valid GKM graph.
"""
function GKMsubgraph_from_edges(gkm::AbstractGKM_graph, edges::Vector{Edge}) :: AbstractGKM_subgraph
  
  vDict = zeros(Int64, 0)

  for e in edges
    @req has_edge(gkm.g, e) "Edge $e not found in GKM graph"
    if !(src(e) in vDict)
      push!(vDict, src(e))
    end
    if !(dst(e) in vDict)
      push!(vDict, dst(e))
    end
  end

  subnv = length(vDict)
  labels = [gkm.labels[vDict[i]] for i in 1:subnv]
  subGKM = gkm_graph(Graph{Undirected}(subnv), labels, gkm.M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}(); checkLabels=false) # use the same character lattice
  
  for e in edges
    sd = indexin([src(e), dst(e)], vDict)
    GKMadd_edge!(subGKM, sd[1], sd[2], gkm.w[e])
  end
  return AbstractGKM_subgraph(gkm, subGKM, vDict)
end

"""
Return true if the gkm subgraph contains the given edge of the supergraph.
"""
function has_edge(gkmSub::AbstractGKM_subgraph, e::Edge)::Bool
  if !(src(e) in gkmSub.vDict) || !(dst(e) in gkmSub.vDict)
    return false
  end
  sd = indexin([src(e), dst(e)], gkmSub.vDict)
  return has_edge(gkmSub.self.g, Edge(sd[1], sd[2]))
end

"""
Return true if the gkm subgraph contains the given vertex of the supergraph.
"""
function has_vertex(gkmSub::AbstractGKM_subgraph, v::Int64)::Bool
  return v in gkmSub.vDict
end

"""
Return true if the gkm subgraph contains the given vertex label of the supergraph.
"""
function has_vertex(gkmSub::AbstractGKM_subgraph, vertexLabel::String)::Bool
  @req vertexLabel in gkmSub.super.labels "Vertex with label $vertexLabel does not exist."
  v::Int64 = indexin([vertexLabel], gkmSub.super.labels)[1]
  return has_vertex(gkmSub, v)
end

function vertexToSupgraph(gkmSub::AbstractGKM_subgraph, v::Int64)::Int64
  @req v in 1:n_vertices(gkmSub.self.g) "Vertex $v not in subgraph."
  return gkmSub.vDict[v]
end

function edgeToSupergraph(gkmSub::AbstractGKM_subgraph, e::Edge)::Edge
  @req has_edge(gkmSub.self.g, e) "Edge $e not contained in subgraph."
  imE = Edge(vertexToSupgraph(gkmSub, src(e)), vertexToSupgraph(gkmSub, dst(e)))
  return imE
end

"""
Return true if the connection map sends edge pairs contained in the subgraph to an edge of the subgraph.
"""
function isCompatible(gkmSub::AbstractGKM_subgraph, con::GKM_connection; printDiagnostics::Bool=true)::Bool
  nvsub = n_vertices(gkmSub.self.g)
  for v in 1:nvsub
    for w in 1:nvsub
      e = Edge(v, w)
      if !has_edge(gkmSub.self.g, e)
        continue
      end
      for u in all_neighbors(gkmSub.self.g, v)
        ei = Edge(v, u)
        epi = con.con[(edgeToSupergraph(gkmSub, e), edgeToSupergraph(gkmSub, ei))]
        if !has_edge(gkmSub, epi)
          printDiagnostics && println("Connection sends image of ($e, $ei) in supergraph to outside the subgraph.")
          return false
        end
      end
    end
  end
  return true
end

function Base.show(io::IO, G::AbstractGKM_subgraph)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM subgraph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM subgraph with $(n_vertices(G.self.g)) nodes and valency $(valency(G.self))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::AbstractGKM_subgraph)

  println(io, "GKM subgraph of:")
  show(io, MIME"text/plain"(), G.super)
  println(io, "\nSubgraph:")
  show(io, MIME"text/plain"(), G.self)
end

"""
Return true if the given GKM subgraph is valid. This holds if and only if all of the following hold:
  1. The supergraph and subgraph are both valid GKM GKMsubgraphs of the same character group
  2. The subgraph is mathematically a subgraph of the supergraph
  3. The edge weights of the subgraph match that of the supergraph
  4. The vertex labels of the subgraph and the supergraph match.
"""
function GKM_isValidSubgraph(gkmsub::AbstractGKM_subgraph; printDiagnostics::Bool = true)::Bool
  if !GKM_isValid(gkmsub.super; printDiagnostics)
    printDiagnostics && println("GKM-Supergraph is invalid")
    return false
  elseif !GKM_isValid(gkmsub.self; printDiagnostics)
    printDiagnostics && println("Sub-GKM-graph is invalid as GKM graph")
    return false
  elseif gkmsub.self.M != gkmsub.super.M
    printDiagnostics && println("GKM parent and subgraph don't have the same character group")
    return false
  end
  
  parentVertices = 1:n_vertices(gkmsub.super.g)
  for v in gkmsub.vDict
    if !(v in parentVertices)
      printDiagnostics && println("Vertex $v not in parent GKM graph")
      return false
    end
  end
  for e in edges(gkmsub.self.g)
    targetEdge = Edge(gkmsub.vDict[src(e)], gkmsub.vDict[dst(e)])
    if !has_edge(gkmsub.super.g, targetEdge)
      printDiagnostics && println{"Edge $e in gets mapped to non-existent edge $targetEdge in parent GKM graph"}
      return false
    elseif gkmsub.self.w[e] != gkmsub.super.w[targetEdge]
      printDiagnostics && println("Weights of $e and its image $targetEdge in the parent GKM graph don't match")
      return false
    end
  end
  
  for v in 1:n_vertices(gkmsub.self.g)
    if gkmsub.self.labels[v] != gkmsub.super.labels[gkmsub.vDict[v]]
      printDiagnostics && println("Label of vertex $v disagrees in subgraph and supergraph.")
      return false
    end
  end
  return true
end