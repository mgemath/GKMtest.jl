using Test

#some working examples
M = free_module(ZZ, 2);
g = Graph{Undirected}(3);
add_edge!(g, 1, 2);
add_edge!(g, 2, 3);
de = Dict(edges(g) .=> [zero(M) for _ in n_edges(g)]);
labs = Symbol.(["a", "b", "c"]);
G = gkm_graph(g, labs, M, de)
[G, G]

de2 = Dict(edges(g) .=> [gens(M)[1] for _ in n_edges(g)]);
G2 = gkm_graph(g, labs, M, de2)

GKMproj_space(4)
GKMproj_space(4, label = "a")
@test true