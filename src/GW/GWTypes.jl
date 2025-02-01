struct GW_decorated_tree
    gkm::AbstractGKM_graph
    tree::Graph
    vDict::Dict{Int, Int} # map vertices of tree to vertices of gkm.g
    edgeMult::Dict{Edge, Int} # each edge of tree has a non-negative multiplicity
    marks::Vector{Int} # vector of marked vertices of the tree 
    R
    
    function GW_decorated_tree(
      gkm::AbstractGKM_graph,
      tree::Graph{Undirected},
      vDict::Dict{Int, Int},
      edgeMult::Dict{Edge, Int},
      marks::Vector{Int},
      R
    )
      return new(gkm, tree, vDict, edgeMult, marks, R)
    end
  end