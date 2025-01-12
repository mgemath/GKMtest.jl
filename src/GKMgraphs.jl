function gkm_graph(
  g::Graph,
  labels::Vector{String},
  M::AbstractAlgebra.Generic.FreeModule{ZZRingElem}, # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}; 
  check::Bool=true
)
  # construct the GKM_graph
  if check
    @req length(labels) == n_vertices(g) "The number of labels does not match the number of fixed points"
    @req all(e -> parent(w[e]) === M, edges(g)) "Character group mismatch"
    @req Set(edges(g)) == keys(w) "The axial function is not well defined"
    @req all(v -> length(all_neighbors(g, 1)) == length(all_neighbors(g, v)), 2:n_vertices(g)) "The valency is not the same for all vertices"
  end
  for e in edges(g)
    w[reverse(e)] = -w[e]
  end
  return AbstractGKM_graph(g, labels, M, w)
end

function valency(G::AbstractGKM_graph)
  return length(all_neighbors(G.g, 1))
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
  # @req valency(G) > 1 "valency is too low"

  # for v in 1:n_vertices(G.g)
  #   for (a, b) in Iterators.product(all_neighbors(G.g, v), all_neighbors(G.g, v))
  #     (a >= b) && continue

  #     if rank(matrix([G.w[Edge(v, a)]; G.w[Edge(v, b)]])) < 2
  #       return false
  #     end

  #   end
  # end

  # return true

end

function is3_indep(G::AbstractGKM_graph)
  return _indep(G, 3)
  # @req valency(G) > 2 "valency is too low"
  
  # for v in 1:n_vertices(G.g)
  #   for (a, b, c) in Iterators.product(all_neighbors(G.g, v), all_neighbors(G.g, v), all_neighbors(G.g, v))
  #     (a >= b || b >= c) && continue

  #     if rank(matrix([G.w[Edge(v, a)]; G.w[Edge(v, b)]; G.w[Edge(v, c)]])) < 3
  #       return false
  #     end

  #   end
  # end

  # return true
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


function GKMadd_edge!(G::AbstractGKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem})

  @req (s in G.labels) "Source not found"
  @req (d in G.labels) "Destination not found"
  @req parent(weight) === G.M "The group of characters is not correct"

  sd = indexin([s, d], G.labels)

  Oscar.add_edge!(G.g, sd[1], sd[2])
  G.w[Edge(sd[1], sd[2])] = weight
  G.w[Edge(sd[2], sd[1])] = -weight

end

function empty_gkm_graph(n::Int, val::Int, labels::Vector{String})

  return gkm_graph(Graph{Undirected}(n), labels, free_module(ZZ, val), Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}())
end