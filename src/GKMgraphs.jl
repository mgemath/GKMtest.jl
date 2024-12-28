function gkm_graph(
  g::Graph{Directed},
  labels::Vector,
  α::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{C}};
  check::Bool=true
) where {C<:RingElement}
  # construct the GKM_graph
  if check
    @req length(labels) == n_vertices(g) "The number of labels does not match the number of fixed points"
    # check 
  end
  return AbstractGKM_graph{C}(g, labels, α; check)
end

function valency(G::AbstractGKM_graph)
  return length(all_neighbors(G.g, 1))
end


function Base.show(io::IO, G::AbstractGKM_graph{C}) where {C<:RingElement}

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM graph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM graph with $(n_vertices(G.g)) nodes and valency $(valency(G))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::AbstractGKM_graph{C}) where {C<:RingElement}

  print(io, "GKM graph with $(n_vertices(G.g)) nodes, valency $(valency(G)) and axial function:")
  for e in edges(G.g)
    print(io, "\n$(G.labels[src(e)]) -> $(G.labels[dst(e)]) => $(G.α[e])")
  end
end


function equivariant_cohomology_ring(G::AbstractGKM_graph{C}) where {C<:RingElement}
    
end