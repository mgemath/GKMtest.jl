using Test

M = free_module(QQ, 2);
g = Graph{Directed}(3);
add_edge!(g, 1, 2);
add_edge!(g, 2, 3);
de = Dict(edges(g) .=> [zero(M) for _ in n_edges(g)]);
labs = ["a", 2, "c"];
G = gkm_graph(g, labs, de)
[G, G]
gkm_graph(g, [1,2], de)


@test true