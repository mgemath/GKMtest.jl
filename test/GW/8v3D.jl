# Implementing the GKM graph from [GKZ20, after Cor 2.13]
# It has a 2-torus acting and is 3-valent, so not 3-independent and admits multiple connections.

g = Graph{Undirected}(8)
labels = ["v_$i" for i in 1:8]
M = free_module(ZZ, 2)
(t1, t2) = gens(M)
w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
G = gkm_graph(g, labels, M, w)

weights = [t1, t2, -t1, -t2]
for i in 1:8
  GKMadd_edge!(G, i, i%8 + 1, weights[(i-1) % 4 + 1])
end
GKMadd_edge!(G, 1, 5, t1+t2)
GKMadd_edge!(G, 2, 6, -t1+t2)
GKMadd_edge!(G, 3, 7, -t1-t2)
GKMadd_edge!(G, 4, 8, -t1-t2)

H2 = GKM_second_homology(G)

#TODO: manually put a connection in place, so that we can calculate QH!

S = QH_structure_constants(G)
pts = [QH_class(G, pointClass(i, G)) for i in 1:8]