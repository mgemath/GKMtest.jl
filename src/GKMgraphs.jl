function gkm_graph(
  g::Graph{Directed},
  α::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{C}};
  check::Bool=true
) where {C<:RingElement}
  # construct the GKM_graph
  if check
    # check 
  end
  return AbstractGKM_graph{C}(g, α; check)
end

function valency(G::AbstractGKM_graph)
  return length(all_neighbors(G.g, 1))
end


function Base.show(io::IO, G::AbstractGKM_graph{C}) where {C<:RingElement}

  print(io, "GKM graph with ", n_vertices(G.g), " nodes, valency ", valency(G), " and axial function: \n", G.α)
end

function equivariant_cohomology_ring(G::AbstractGKM_graph{C}) where {C<:RingElement}
    
end