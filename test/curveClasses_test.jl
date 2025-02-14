# Projective space example:
P3 = GKMproj_space(3)
H2 = GKM_second_homology(P3)
beta = 5*gens(H2.H2)[1]
edgeList = [Edge(1, 2), Edge(2, 3), Edge(2,1)]

ms = GKMtest._multiplicities(H2, edgeList, beta)

println("Multiplicities for P3 with beta=$beta and edgeList: $edgeList")
println("Dual cone ray sum: $(H2.dualConeRaySum)")
println("Max number of edges: $(GKMtest._max_n_edges(H2, beta))")
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
edgeList = Edge[]
push!(edgeList, edgeFromLabels(X, "x_0,x_0", "x_1,x_0")) # the first three edges represent c1
push!(edgeList, edgeFromLabels(X, "x_0,x_0", "x_2,x_0"))
push!(edgeList, edgeFromLabels(X, "x_1,x_0", "x_2,x_0"))
push!(edgeList, edgeFromLabels(X, "x_0,x_0", "x_0,x_1")) # this edge represents c2

ms = GKMtest._multiplicities(H2, edgeList, beta)

println("Multiplicities for P2xP2 with beta=$beta and edgeList: $edgeList")
println("Dual cone ray sum: $(H2.dualConeRaySum)")
println("Max number of edges: $(GKMtest._max_n_edges(H2, beta))")
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
edgeList = Edge[]
push!(edgeList, edgeFromLabels(F, "412", "142")) # this edge represents a1
push!(edgeList, edgeFromLabels(F, "412", "421")) # this edge represents a2
push!(edgeList, edgeFromLabels(F, "412", "413")) # this edge represents a3
push!(edgeList, edgeFromLabels(F, "234", "231")) # this edge also represents a3

ms = GKMtest._multiplicities(H2, edgeList, beta)

println("Multiplicities for Flag(C^4) with beta=$beta and edgeList: $edgeList")
println("Dual cone ray sum: $(H2.dualConeRaySum)")
println("Max number of edges: $(GKMtest._max_n_edges(H2, beta))")
for m in ms
  println(m)
end

# P3 blown up somewhere
P3 = projective_space(NormalToricVariety,3);#thespaceP^3.
X = gkm_graph(domain(blow_up(P3,ideal(gens(cox_ring(P3))[1:3]))));
H2 = GKM_second_homology(X);
R = equivariant_cohomology_ring(X);

println("Classes of edges of blown-up P3")
for e in edges(X.g)
    println("$(X.labels[src(e)]) -> $(X.labels[dst(e)]) => ", edgeCurveClass(H2, e))
end

H=edgeCurveClass(H2, edgeFromLabels(X, "5", "4")); #(0,1)
# E=edgeCurveClass(H2, edgeFromLabels(X, "1", "4")) #(1,0)
E=edgeCurveClass(H2, edgeFromLabels(X, "1", "2")) #(-1,1)
d=2;#thiscanbeanynonnegativeinteger
e=-2;#thiscanbeanynonpositiveinteger
beta=d*H+e*E;
# beta = 3*E

edgeList = Edge[]
push!(edgeList, edgeFromLabels(X, "1", "4")) # this edge represents (1,0)
push!(edgeList, edgeFromLabels(X, "1", "4")) # this edge represents (1,0)
push!(edgeList, edgeFromLabels(X, "1", "2")) # this edge represents (-1,1)

ms = GKMtest._multiplicities(H2, edgeList, beta)

println("Multiplicities for Bl(P3) with beta=$beta and edgeList: $edgeList")
println("Dual cone ray sum: $(H2.dualConeRaySum)")
println("Max number of edges: $(GKMtest._max_n_edges(H2, beta))")
for m in ms
  println(m)
end
