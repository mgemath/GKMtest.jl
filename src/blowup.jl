"""
Return (GKM graph of blowup, GKM graph of exceptional divisor, connection on blown-up graph)
from (GKM graph, GKM subgraph, connection on supergraph), where both are encoded as AbstractGKM_subgraph.
This follows [Guillemin--Zara, section 2.2.1]

Warning: This will build an Undirected graph. Behaviour with directed graphs as input is not tested.
"""
function blowupGKM(gkmSub::AbstractGKM_subgraph, con::GKM_connection)::Tuple{AbstractGKM_subgraph, GKM_connection}
  
  @req GKM_isValidSubgraph(gkmSub) "invalid graph/subgraph pair"
  @req GKM_isValidConnection(con) "invalid connection"
  @req isCompatible(gkmSub, con) "connection incompatible with subgraph"

  nvSub = n_vertices(gkmSub.self.g)
  nv = n_vertices(gkmSub.super.g)
  vDict = gkmSub.vDict
  M = gkmSub.super.M
  d = valency(gkmSub.self)
  n = valency(gkmSub.super)
  c = n - d # codimension

  if d == n
    return (GKMsubgraph_from_vertices(gkmSub.super, Array(1:nv)), con)
  end
  
  externalVertices = Int64[]
  for i in 1:nv
    if !has_vertex(gkmSub, i)
      push!(externalVertices, i)
    end
  end

  # to each vertex of the subgraph, associate its neighbors in the supergraph s.t. the connecting edge is not in the subgraph
  normalNeighbors = Vector{Int64}[]
  for i in 1:nvSub
    tmp = Int64[]
    for j in all_neighbors(gkmSub.super.g, vDict[i])
      if !has_edge(gkmSub, Edge(vDict[i], j))
        push!(tmp, j)
      end
    end
    push!(normalNeighbors, tmp)
  end

  nvBlowup = (c * nvSub) + (nv - nvSub)
  labels = String[]

  for i in 1:nvSub
    for j in normalNeighbors[i]
      push!(labels, "[" * gkmSub.self.labels[i] * ">" * gkmSub.super.labels[j] * "]")
    end
  end

  for v in externalVertices
    push!(labels, gkmSub.super.labels[v])
  end

  gkmBlowup = gkm_graph(Graph{Undirected}(nvBlowup), labels, M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}(); checkLabels=false)
  exceptionalEdges = Edge[]
  blowupCon = Dict{Tuple{Edge, Edge}, Edge}()

  # build complete graph for each blown up vertex
  for v in 1:nvSub
    for i in normalNeighbors[v], j in normalNeighbors[v]
      if i >= j
        continue
      end
      weight = gkmSub.super.w[Edge(vDict[v], j)] - gkmSub.super.w[Edge(vDict[v], i)]
      newEdge = Edge(intVindex(v, c, i, normalNeighbors), intVindex(v, c, j, normalNeighbors))
      GKMadd_edge!(gkmBlowup, src(newEdge), dst(newEdge), weight)
      push!(exceptionalEdges, newEdge)

      # println("connection for new edge $newEdge:")
      # build connection for newEdge & its reverse
      for k in normalNeighbors[v]
        if k == i
          eiDst = extFlagToIndex(gkmSub, Edge(i, vDict[v]), normalNeighbors, externalVertices, nvBlowup, c)
          epiDst = extFlagToIndex(gkmSub, Edge(j, vDict[v]), normalNeighbors, externalVertices, nvBlowup, c)
          ei = Edge(src(newEdge), eiDst)
          epi = Edge(dst(newEdge), epiDst)
          blowupCon[(newEdge, ei)] = epi
          blowupCon[(reverse(newEdge), epi)] = ei
          # println("($newEdge, $ei) -> $epi")
        elseif k == j
          blowupCon[(newEdge, newEdge)] = reverse(newEdge)
          blowupCon[(reverse(newEdge), reverse(newEdge))] = newEdge
          # println("($newEdge, $newEdge) -> $(reverse(newEdge))")
        else
          vertKInd = intVindex(v, c, k, normalNeighbors)
          ei = Edge(src(newEdge), vertKInd)
          epi = Edge(dst(newEdge), vertKInd)
          blowupCon[(newEdge, ei)] = epi
          blowupCon[(reverse(newEdge), epi)] = ei
          # println("($newEdge, $ei) -> $epi")
        end
      end
      
      # build connection for (newEdge, ei) where ei connects to an internal neighbor of v.
      for n in all_neighbors(gkmSub.self.g, v)
        vn = Edge(vDict[v], vDict[n])
        vi = Edge(vDict[v], i)
        vj = Edge(vDict[v], j)
        na = con.con[(vn, vi)]
        nb = con.con[(vn, vj)]
        a = dst(na)
        b = dst(nb)
        naInd = intVindex(n, c, a, normalNeighbors)
        nbInd = intVindex(n, c, b, normalNeighbors)
        ei = Edge(src(newEdge), naInd)
        epi = Edge(dst(newEdge), nbInd)
        blowupCon[(newEdge, ei)] = epi
        blowupCon[(reverse(newEdge), epi)] = ei
        #println("($newEdge, $ei) -> $epi")
      end
    end
  end

  # println("Connection from old edges:")

  # need this data to translate connection to new edges coming from original edges
  flagBij = flagBijections(gkmSub, con, normalNeighbors, c, externalVertices, nvBlowup)

  # add edges coming from original edges
  for e in edges(gkmSub.super.g)
    s = src(e)
    d = dst(e)
    
    if !has_edge(gkmSub, e)

      sInd = extFlagToIndex(gkmSub, e, normalNeighbors, externalVertices, nvBlowup, c)
      dInd = extFlagToIndex(gkmSub, reverse(e), normalNeighbors, externalVertices, nvBlowup, c)
      w = gkmSub.super.w[e]
      GKMadd_edge!(gkmBlowup, sInd, dInd, w)
      addToConnection!(blowupCon, con, e, Edge(sInd, dInd), flagBij)
    else

      sSubInd = indexin(s, vDict)[1]

      for n in normalNeighbors[sSubInd]

        ei = Edge(s, n)
        epi = con.con[(e, ei)]
        sIndi = extFlagToIndex(gkmSub, ei, normalNeighbors, externalVertices, nvBlowup, c)
        dIndi = extFlagToIndex(gkmSub, epi, normalNeighbors, externalVertices, nvBlowup, c)
        w = gkmSub.super.w[e]

        GKMadd_edge!(gkmBlowup, sIndi, dIndi, w)
        addToConnection!(blowupCon, con, e, Edge(sIndi, dIndi), flagBij)
        push!(exceptionalEdges, Edge(sIndi, dIndi))
      end
    end
  end

  gkmSubgraphBlowup = GKMsubgraph_from_edges(gkmBlowup, exceptionalEdges)
  # Base.show(stdout, MIME"text/plain"(), gkmSubgraphBlowup)
  #println("Resulting connection dict:")
  for k in keys(blowupCon)
    #println("$k -> $(blowupCon[k])")
  end
  fullBlowupConnection = build_GKM_connection(gkmBlowup, blowupCon)

  return (gkmSubgraphBlowup, fullBlowupConnection)

end

function addToConnection!(blowupCon::Dict{Tuple{Edge, Edge}, Edge}, con::GKM_connection, oldE::Edge, newE::Edge, flagBij::Array{Dict{Int64, Int64}})
  newS = src(newE)
  newD = dst(newE)
  for (e,ei) in keys(con.con)
    e != oldE && continue
    epi = con.con[(e, ei)]
    newEi = Edge(newS, flagBij[newS][dst(ei)])
    newEpi = Edge(newD, flagBij[newD][dst(epi)])
    blowupCon[(newE, newEi)] = newEpi
    blowupCon[(reverse(newE), newEpi)] = newEi
    # println("($newE, $newEi) -> $newEpi")
  end
end

"""
For each vertex of the blowup graph, this returns a dictionary from neighbors in the original graph to neighbors in the new graph.
"""
function flagBijections(gkmSub::AbstractGKM_subgraph, con::GKM_connection, normalNeighbors::Vector{Vector{Int64}}, c::Int64, externalVertices::Vector{Int64}, nvBlowup::Int64)

  res = Array{Dict{Int64, Int64}}(undef, nvBlowup)
  vDict = gkmSub.vDict

  # calculate external vertices' flags
  for v in externalVertices
    resV = Dict{Int64, Int64}()
    for n in all_neighbors(gkmSub.super.g, v)
      nInd = extFlagToIndex(gkmSub, Edge(n, v), normalNeighbors, externalVertices, nvBlowup, c)
      resV[n] = nInd
    end
    vInd = extVindex(v, nvBlowup, externalVertices)
    res[vInd] = resV
  end

  # calculate internal vertices' flags
  for v in 1:n_vertices(gkmSub.self.g)
    for n in normalNeighbors[v]
      resVn = Dict{Int64, Int64}()
      for w in normalNeighbors[v]
        if n == w
          resVn[w] = extFlagToIndex(gkmSub, Edge(n, vDict[v]), normalNeighbors, externalVertices, nvBlowup, c)
        else
          resVn[w] = intVindex(v, c, w, normalNeighbors)
        end
      end
      for u in all_neighbors(gkmSub.self.g, v)
        e = Edge(vDict[v], vDict[u])
        ei = Edge(vDict[v], n)
        epi = con.con[(e, ei)]
        m = dst(epi)
        resVn[vDict[u]] = intVindex(u, c, m, normalNeighbors)
      end
      vnInd = intVindex(v, c, n, normalNeighbors)
      res[vnInd] = resVn
    end
  end

  return res
end

function intVindex(v::Int64, c::Int64, i::Int64, normalNeighbors::Vector{Vector{Int64}})
  return (v-1)*c + indexin(i, normalNeighbors[v])[1]
end

function extVindex(v::Int64, nvBlowup::Int64, externalVertices::Vector{Int64})
  i::Int64 =  indexin(v, externalVertices)[1]
  return nvBlowup - length(externalVertices) + i
end

"""
For an external flag at src(e), return the index of the source of that flag in the bowup.
"""
function extFlagToIndex(gkmSub::AbstractGKM_subgraph, e::Edge, normalNeighbors::Vector{Vector{Int64}}, externalVertices::Vector{Int64}, nvBlowup::Int64, c::Int64)::Int64
  s = src(e)
  d = dst(e)
  if has_vertex(gkmSub, s)
    sInSub = indexin(s, gkmSub.vDict)[1]
    sInd = intVindex(sInSub, c, d, normalNeighbors)
    return sInd
  else
    return extVindex(s, nvBlowup, externalVertices)
  end
end