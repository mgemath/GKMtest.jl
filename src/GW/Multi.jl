#TODO implement an Iterator of multiplicities

function all_ones_multip(G::AbstractGKM_graph, tree::Graph)

    return [Dict{Edge, Int}(edges(tree) .=> [1 for _ in 1:n_edges(tree)])]
end