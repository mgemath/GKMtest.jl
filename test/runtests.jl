using Test

M = free_module(QQ, 2);
g = Graph{Directed}(4);
add_edge!(g, 1, 2);
add_edge!(g, 2, 3);
de = Dict(edges(g) .=> [zero(M) for _ in n_edges(g)])
gkm_graph(g, de)

@test true