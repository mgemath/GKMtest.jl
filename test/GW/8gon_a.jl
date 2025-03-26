# Implementing the left hand GKM graph from [GKZ22, Example 2.44]
# As unsigned GKM graph, it is realized by the equivariant connected sum of three copies of S^2xS^2

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

H2 = GKM_second_homology(G)

S = QH_structure_constants(G)
pts = [QH_class(G, pointClass(i, G)) for i in 1:8]
edg = [QH_class(G, PDClass(GKMsubgraph_from_vertices(G, [i, i%8 + 1]))) for i in 1:8]

# julia> QH_is_commutative(G)
#true
#
#julia> QH_is_associative(G)
#Not associative for: (1, 1, 2)
#false

# (Q:) is this just because we didn't check high enough q coefficients?
# (A:) I guess calculating higher chern number stuff won't cancel out the non-associativity in the chern class 2 part!
# So this is actually NOT associative!