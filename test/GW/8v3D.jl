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

conDict1 = Dict{Tuple{Edge, Edge}, Edge}(
    (Edge(1, 2), Edge(1, 5)) => Edge(2, 6),
    (Edge(1, 2), Edge(1, 8)) => Edge(2, 3),
    (Edge(2, 3), Edge(2, 1)) => Edge(3, 4),
    (Edge(2, 3), Edge(2, 6)) => Edge(3, 7),
    (Edge(3, 4), Edge(3, 7)) => Edge(4, 8),
    (Edge(3, 4), Edge(3, 2)) => Edge(4, 5),
    (Edge(4, 5), Edge(4, 3)) => Edge(5, 6),
    (Edge(4, 5), Edge(4, 8)) => Edge(5, 1),
    (Edge(5, 6), Edge(5, 1)) => Edge(6, 2),
    (Edge(5, 6), Edge(5, 4)) => Edge(6, 7),
    (Edge(6, 7), Edge(6, 2)) => Edge(7, 3),
    (Edge(6, 7), Edge(6, 5)) => Edge(7, 8),
    (Edge(7, 8), Edge(7, 3)) => Edge(8, 4),
    (Edge(7, 8), Edge(7, 6)) => Edge(8, 1),
    (Edge(8, 1), Edge(8, 7)) => Edge(1, 2),
    (Edge(8, 1), Edge(8, 4)) => Edge(1, 5),
    # Now the diagonal edges:
    (Edge(1, 5), Edge(1, 2)) => Edge(5, 6),
    (Edge(1, 5), Edge(1, 8)) => Edge(5, 4),
    (Edge(2, 6), Edge(2, 1)) => Edge(6, 5),
    (Edge(2, 6), Edge(2, 3)) => Edge(6, 7),
    (Edge(3, 7), Edge(3, 2)) => Edge(7, 6),
    (Edge(3, 7), Edge(3, 4)) => Edge(7, 8),
    (Edge(4, 8), Edge(4, 3)) => Edge(8, 7),
    (Edge(4, 8), Edge(4, 5)) => Edge(8, 1)
)
for e in edges(G.g)
  conDict1[(e, e)] = reverse(e)
end
for (e, ei) in keys(conDict1)
  epi = conDict1[(e, ei)]
  conDict1[reverse(e), epi] = ei
end
con1 = build_GKM_connection(G, conDict1)
set_GKM_connection!(G, con1)

# Generators (0, -1, 1), (0, 1, 0), (1, -1, 0) of Chern classes 0, 2, 0, resp.
H2gens = [edgeCurveClass(G, e) for e in [Edge(5, 6), Edge(3, 7), Edge(6, 7)]]

# This takes long, but one can just abort whenever, and all betas that have been calculated will be stored in G.QH_structure_constants.
for i in 0:3, j in 0:3, k in 0:3
  b = i*H2gens[1] + j*H2gens[2] + k*H2gens[3]
  QH_structure_constants(G, b)
end

# Result:
# We get polynomial QH structure constants in all curve classes computed so far.
# The nonzero ones calculated so far are:
#(0, -5, 5)
#(0, 1, 0)
#(1, -1, 0)
#(1, 0, 0)
#(1, 1, 0)
#(3, -3, 0)
#(0, -4, 4)
#(0, 0, 0)
#(2, -2, 0)
#
# Furthermore, the coefficients for beta=(0, -j, j) are constant for j in {2, 3, 4, 5}, so probably we get some q^beta/(1-q^beta) terms!
# The same holds for beta=(j, -j, 0) for j in {1, 2, 3}
# Note that the other classes can not appear infinitely often (if the structure constants are polynomial) because they have positive chern number.