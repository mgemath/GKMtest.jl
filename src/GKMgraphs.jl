function gkm_graph(
  g::Graph,
  labels::Vector{Symbol},
  M::AbstractAlgebra.Generic.FreeModule{ZZRingElem}, # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}; 
  check::Bool=true
)
  # construct the GKM_graph
  if check
    @req length(labels) == n_vertices(g) "The number of labels does not match the number of fixed points"
    @req all(e -> parent(w[e]) === M, edges(g)) "Character group mismatch"
    @req Set(edges(g)) == keys(w) "The axial function is not well defined"
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


function equivariant_cohomology_ring(G::AbstractGKM_graph)
    
end

function rank_torus(G::AbstractGKM_graph)
  return rank(G.M)
end

function connection(G::AbstractGKM_graph, e::Edge)

  # return a dictionary edges(g) -> (Edge, ZZ) where for each edge e_i containing src(e), connection[e_i] = (e', a_i) such that 
  # w[e', dst(e)] = w[e_i, src(e)] - a_i * w[e, src(e)]
  return   
end

function GKMproj_space(dim::Int; label::String = "x_")
  g = complete_graph(dim+1)
  M = free_module(ZZ, dim+1)
  w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for e in edges(g)
    w[e] = gens(M)[src(e)]-gens(M)[dst(e)]
  end
  labels = [Symbol(label*"$i") for i in 0:dim]
  return gkm_graph(g, labels, M, w)
end