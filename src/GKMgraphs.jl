"""
Create a GKM graph from the given data.
Warnings:
 1. Do not change the number of verices after this.
 2. After you have added all edges using GKMadd_edge!(), use GKM_initialize!() to calculate
    the GKM connection (if it is unique) and the curve classes.
"""
function gkm_graph(
  g::Graph,
  labels::Vector{String},
  M::AbstractAlgebra.Generic.FreeModule{R}, # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}}; 
  check::Bool=true,
  checkLabels::Bool=true
) where R <: GKM_weight_type
  # construct the GKM_graph
  if check
    @req n_vertices(g) >= 1 "GKM graph needs at least one vertex"
    @req length(labels) == n_vertices(g) "The number of labels does not match the number of fixed points"
    @req all(e -> parent(w[e]) === M, edges(g)) "Character group mismatch"
    @req Set(edges(g)) == keys(w) "The axial function is not well defined"
    @req all(v -> length(all_neighbors(g, 1)) == length(all_neighbors(g, v)), 2:n_vertices(g)) "The valency is not the same for all vertices"
    @req length(unique(labels)) == length(labels) "Labels must be unique"
  end
  if checkLabels
    # reserve characters <,[,] for vertex labels of blowups
    @req all(v -> !contains(labels[v], ">") && !contains(labels[v], "[") && !contains(labels[v], "]"), 1:n_vertices(g)) "Characters >,[,] are forbidden for vertex labels"
  end
  for e in edges(g)
    w[reverse(e)] = -w[e]
  end

  GW_structure_consts = Dict{CurveClass_type, Array{Any, 3}}()

  gkm = AbstractGKM_graph(g, labels, M, w, nothing, nothing, nothing, GW_structure_consts, false)
  gkm.equivariantCohomology = _equivariant_cohomology_ring(gkm)

  return gkm
end

"""
Call this as soon as all edges have been added to the GKM graph to calculate the GKM connection
(if unique) and the curve classes of the gkm graph.
This will set the fields gkm.connection and gkm.curveClasses that are initially nothing.
If you don't call this, these fields will be initialized later if possible, which might take some time
(especially for curveClasses).
If any of those fields are already set, this will not overwrite them.
"""
function GKM_initialize!(gkm::AbstractGKM_graph; connection::Bool=true, curveClasses::Bool=true)

  if connection
    get_GKM_connection(gkm)
  end
  if curveClasses
    GKM_second_homology(gkm)
  end
end

function valency(G::AbstractGKM_graph)
  return length(all_neighbors(G.g, 1))
end

function _common_weight_denominator(G::AbstractGKM_graph)::ZZRingElem
  if G.weightType <: ZZRingElem
    return ZZ(1)
  elseif G.weightType <: QQFieldElem
    res::ZZRingElem = ZZ(1)
    rk = rank_torus(G)
    for e in edges(G.g)
      for i in 1:rk
        res = lcm(res, denominator(G.w[e][i]))
      end
    end
    return res
  end
  @req false "Ony ZZRingElem and QQFieldElem are supported as weight types so far."
end

function Base.show(io::IO, G::AbstractGKM_graph)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM graph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM graph with $(n_vertices(G.g)) nodes and valency $(valency(G))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::AbstractGKM_graph)

  print(io, "GKM graph with $(n_vertices(G.g)) nodes, valency $(valency(G)) and axial function:")
  for e in edges(G.g)
    print(io, "\n$(G.labels[src(e)]) -> $(G.labels[dst(e)]) => $(G.w[e])")
  end
end

function rank_torus(G::AbstractGKM_graph)
  return rank(G.M)
end


function GKMproj_space(dim::Int; label::String = "x_")
  g = complete_graph(dim+1)
  M = free_module(ZZ, dim+1)
  w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for e in edges(g)
    w[e] = gens(M)[src(e)]-gens(M)[dst(e)]
  end
  labels = [label*"$i" for i in 0:dim]
  return gkm_graph(g, labels, M, w)
end

function is2_indep(G::AbstractGKM_graph)
  return _indep(G, 2)
end

function is3_indep(G::AbstractGKM_graph)
  return _indep(G, 3)
end

function _indep(G::AbstractGKM_graph, k::Int64)
  
  @req valency(G) >= k "valency is too low"

  for v in 1:n_vertices(G.g)
    for tup in Iterators.product([all_neighbors(G.g, v) for _ in 1:k]...)
      any(i-> tup[i-1] >= tup[i], 2:k) && continue
      
      if rank(matrix([G.w[Edge(v, tup[i])] for i in 1:k])) < k
        return false
      end

    end
  end
  
  return true
end

"""
Warning: Add all edges immediately after creation.
Any curve class or cohomology functionality should only be used after all edges have been added.
The same holds for GKM_initialize().
"""
function GKMadd_edge!(G::AbstractGKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type

  @req (s in G.labels) "Source label not found"
  @req (d in G.labels) "Destination label not found"

  sd = indexin([s, d], G.labels)

  GKMadd_edge!(G, sd[1], sd[2], weight)

end

function GKMadd_edge!(G::AbstractGKM_graph, s::Int64, d::Int64, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type
  
  @req (s in 1:n_vertices(G.g)) "Source $s not found"
  @req (d in 1:n_vertices(G.g)) "Destination $d not found"
  @req parent(weight) === G.M "The group of characters is not correct"

  Oscar.add_edge!(G.g, s, d)
  G.w[Edge(s, d)] = weight
  G.w[Edge(d, s)] = -weight

end

function empty_gkm_graph(n::Int, val::Int, labels::Vector{String})

  return gkm_graph(Graph{Undirected}(n), labels, free_module(ZZ, val), Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}())
end

"""
Return true if the GKM graph is valid. This means:
  1. Every vertex has the same degree
  2. The weights are defined for every edge and every reverse of every edge
  3. The weights belong to the weight lattice
  4. The weights of an edge and its reverse sum to zero
  5. There are the right number of vertex labels
  6. If the valency is at least two, the weights of the graph are 2-independent.
  7. Vertex labels must be unique
  8. The equivariant cohomology ring has rank = number of vertices of graph
  9. The coefficient ring of the equivariant cohomology ring has number of generators = torus rank.
"""
function GKM_isValid(gkm::AbstractGKM_graph; printDiagnostics::Bool=true)::Bool

  if !all(v -> length(all_neighbors(gkm.g, 1)) == length(all_neighbors(gkm.g, v)), 2:n_vertices(gkm.g))
    printDiagnostics && println("The valency is not the same for all vertices")
    return false
  end

  for e in edges(gkm.g)
    if !haskey(gkm.w, e)
      printDiagnostics && println("Weight of $e is missing.")
      return false
    elseif !haskey(gkm.w, reverse(e))
      printDiagnostics && println("Weight of $(reverse(e)) is missing.")
      return false
    elseif parent(gkm.w[e]) != gkm.M
      printDiagnostics && println("Weight of $e doesn't belong to $(gkm.M).")
      return false
    elseif parent(gkm.w[reverse(e)]) != gkm.M
      printDiagnostics && println("Weight of $(reverse(e)) doesn't belong to $(gkm.M).")
      return false
    elseif !(gkm.w[e] == -gkm.w[reverse(e)])
      printDiagnostics && println("Weights of $e and $(reverse(e)) don't sum to zero.")
      return false
    end
  end

  if length(gkm.labels) != n_vertices(gkm.g)
    printDiagnostics && println("Not the right number of labels")
    return false
  elseif (valency(gkm) > 1 && !is2_indep(gkm))
    printDiagnostics && println("GKM graph is not 2-independent.")
    return false
  end

  if length(unique(gkm.labels)) != length(gkm.labels)
    printDiagnostics && println("Labels are not unique.")
    return false
  end

  if (valency(gkm) > 2 && !is3_indep(gkm))
    printDiagnostics && println("GKM graph is valid but not 3-independent, so connections may not be unique.")
  end

  if length(gens(gkm.equivariantCohomology.coeffRing)) != rank_torus(gkm)
    printDiagnostics && println("Coefficient ring of equivariant cohomology has wrong rank.")
    return false
  end

  if length(gens(gkm.equivariantCohomology.cohomRing)) != n_vertices(gkm.g)
    printDiagnostics && println("Equivariant cohomology ring should have rank = number of vertices.")
    return false
  end

  return true
end

function edgeFromLabels(G::AbstractGKM_graph, s::String, d::String)::Edge
  @req s in G.labels "source $s is not a vertex in $G."
  @req d in G.labels "destination $d is not a vertex in $G."

  sd = indexin([s, d], G.labels)
  return Edge(sd[1], sd[2])
end