# Projective space example:
P3 = GKMproj_space(3)
H2 = GKM_second_homology(P3)
beta = 5*gens(H2.H2)[1]
edgeList = [Edge(1, 2), Edge(2, 3), Edge(2,1)]
chern = 4

ms = GKMtest._multiplicities(H2, edgeList, beta)

println("Multiplicities for P3 with beta=$beta and edgeList: $edgeList")
println("Dual cone ray sum: $(H2.dualConeRaySum)")
println("Max number of edges: $(GKMtest._max_n_edges(H2, beta))")
for m in ms
  println(m)
end
println("Dual cone rays: $(rays(H2.dualCone))")
println("Classes with chern class $chern:")
for b in GKMtest._effectiveClassesWithChernNumber(H2, chern)
  println(b)
end


# Product of P2 and P3 example
P2 = GKMproj_space(2)
P3 = GKMproj_space(4)
X = P2*P3
H2 = GKM_second_homology(X)
c1 = edgeCurveClass(H2, edgeFromLabels(X, "x_0,x_0", "x_1,x_0"))
c2 = edgeCurveClass(H2, edgeFromLabels(X, "x_0,x_0", "x_0,x_1"))
beta = 5 * c1 + 7*c2
chern = 6
edgeList = Edge[]
push!(edgeList, edgeFromLabels(X, "x_0,x_0", "x_1,x_0")) # the first three edges represent c1
push!(edgeList, edgeFromLabels(X, "x_0,x_0", "x_2,x_0"))
push!(edgeList, edgeFromLabels(X, "x_1,x_0", "x_2,x_0"))
push!(edgeList, edgeFromLabels(X, "x_0,x_0", "x_0,x_1")) # this edge represents c2

ms = GKMtest._multiplicities(H2, edgeList, beta)

println("Multiplicities for P2xP3 with beta=$beta and edgeList: $edgeList")
println("Dual cone ray sum: $(H2.dualConeRaySum)")
println("Max number of edges: $(GKMtest._max_n_edges(H2, beta))")
for m in ms
  println(m)
end
println("Dual cone rays: $(rays(H2.dualCone))")
println("Classes with chern class $chern:")
for b in GKMtest._effectiveClassesWithChernNumber(H2, chern)
  println(b)
end

# Full flag variety in C^4
F = flag_gkm_graph([1, 1, 1, 1])
H2 = GKM_second_homology(F)
a1 = edgeCurveClass(H2, edgeFromLabels(F, "123", "213"))
a2 = edgeCurveClass(H2, edgeFromLabels(F, "123", "132"))
a3 = edgeCurveClass(H2, edgeFromLabels(F, "123", "124"))
beta = 5*a1 + 2*a2 + 3*a3
edgeList = Edge[]
chern = 6
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
println("Dual cone rays: $(rays(H2.dualCone))")
println("Classes with chern class $chern:")
for b in GKMtest._effectiveClassesWithChernNumber(H2, chern)
  println(b)
end

# P3 blown up somewhere
P3 = projective_space(NormalToricVariety,3)#thespaceP^3.
X = gkm_graph(domain(blow_up(P3,ideal(gens(cox_ring(P3))[1:3]))))
H2 = GKM_second_homology(X)
chern = 4

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
# this is an example where the dual of the effecetive cone is generated by (1,1) and (0,1).
println("Dual cone rays: $(rays(H2.dualCone))")
println("Are the following classes effective?")
for beta in [H, E, H+E, H-E, H-2*E, -H+E]
  println("$beta: $(GKMtest.isEffectiveCurveClass(H2, beta))")
end
println("Classes with chern class $chern:")
for b in GKMtest._effectiveClassesWithChernNumber(H2, chern)
  println(b)
end