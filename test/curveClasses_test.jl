# Projective space example:
P3 = GKMproj_space(3)
H2 = GKM_second_homology(P3)
beta = 5*gens(H2.H2)[1]
edges = [Edge(1, 2), Edge(2, 3), Edge(2,1)]

ms = GKMtest._multiplicities(H2, edges, beta)

println("Multiplicities for P3 with beta=$beta and edges: $edges")
for m in ms
  println(m)
end

# Product of P2 and P2 example
P2 = GKMproj_space(2)
X = P2*P2
H2 = GKM_second_homology(X)
c1 = edgeCurveClass(H2, edgeFromLabels(X, "x_0,x_0", "x_1,x_0"))
c2 = edgeCurveClass(H2, edgeFromLabels(X, "x_0,x_0", "x_0,x_1"))
beta = 5 * c1 + 7*c2
edges = Edge[]
push!(edges, edgeFromLabels(X, "x_0,x_0", "x_1,x_0")) # the first three edges represent c1
push!(edges, edgeFromLabels(X, "x_0,x_0", "x_2,x_0"))
push!(edges, edgeFromLabels(X, "x_1,x_0", "x_2,x_0"))
push!(edges, edgeFromLabels(X, "x_0,x_0", "x_0,x_1")) # this edge represents c2

ms = GKMtest._multiplicities(H2, edges, beta)

println("Multiplicities for P2xP2 with beta=$beta and edges: $edges")
for m in ms
  println(m)
end

# Full flag variety in C^4
F = flag_gkm_graph([1, 1, 1, 1])
H2 = GKM_second_homology(F)
a1 = edgeCurveClass(H2, edgeFromLabels(F, "123", "213"))
a2 = edgeCurveClass(H2, edgeFromLabels(F, "123", "132"))
a3 = edgeCurveClass(H2, edgeFromLabels(F, "123", "124"))
beta = 5*a1 + 2*a2 + 3*a3
edges = Edge[]
push!(edges, edgeFromLabels(F, "412", "142")) # this edge represents a1
push!(edges, edgeFromLabels(F, "412", "421")) # this edge represents a2
push!(edges, edgeFromLabels(F, "412", "413")) # this edge represents a3
push!(edges, edgeFromLabels(F, "234", "231")) # this edge also represents a3

ms = GKMtest._multiplicities(H2, edges, beta)

println("Multiplicities for Flag(C^4) with beta=$beta and edges: $edges")
for m in ms
  println(m)
end