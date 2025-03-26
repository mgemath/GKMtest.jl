# Implementing the right hand GKM graph from [GKZ22, Example 2.44]
# It does not admit a connection (or even an almost complex structure, I think), and cannot be realized by a Hamiltonian T-2 action.
# As unsigned GKM graph, it is realized by the equivariant connected sum of three copies of S^2xS^2

g = Graph{Undirected}(8)
labels = ["v_$i" for i in 1:8]
M = free_module(ZZ, 2)
(t1, t2) = gens(M)
w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
G = gkm_graph(g, labels, M, w)

weights = [t1, t2, t1, -t2, -t1, t2, t1, -t2]
for i in 1:8
  GKMadd_edge!(G, i, i%8 + 1, weights[i])
end

H2 = GKM_second_homology(G) # crashes because c1 is not a GKM class. That's ok.

S = QH_structure_constants(G)