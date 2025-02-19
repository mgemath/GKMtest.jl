struct GW_decorated_tree
    gkm::AbstractGKM_graph
    tree::Graph
    vDict::Union{Vector{Int}, Tuple{Vararg{Int}}} # map vertices of tree to vertices of gkm.g
    edgeMult::Dict{Edge, Int} # each edge of tree has a non-negative multiplicity
    marks::Vector{Int} # vector of marked vertices of the tree 
    
    function GW_decorated_tree(
      gkm::AbstractGKM_graph,
      tree::Graph{Undirected},
      vDict::Union{Vector{Int}, Tuple{Vararg{Int}}},
      edgeMult::Dict{Edge, Int},
      marks::Vector{Int}
    )
      return new(gkm, tree, vDict, edgeMult, marks)
    end
  end