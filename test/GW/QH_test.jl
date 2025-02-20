G = GKMproj_space(3)
#G = blowupGKM(GKMsubgraph_from_vertices(G, [1, 2])).super
beta = edgeCurveClass(G, Edge(1, 2))
class1 = pointClass(1, G)
class2 = pointClass(2, G)
#class3 = pointClass(3, G)
t1, t2 = gens(G.equivariantCohomology.coeffRing)

a = quantumProduct(G, beta, class1, class2)
b = quantumProduct(G, beta, class1, class2; useStructureConstants=false)
#println(a)
#println(b)
println(a==b)
a = quantumProduct(G, beta, class1*class1, one(G.equivariantCohomology))
b = quantumProduct(G, beta, class1*class1, one(G.equivariantCohomology); useStructureConstants=false)
println(a==b)

S = GKMtest.QH_structure_constants(G)
Sfactored = Dict{GKMtest.CurveClass_type, Any}()
for beta in keys(S)
  println("Structure constants for beta=$beta:")
  Sfactored[beta] = [S[beta][k] == 0 ? 0 : factor(numerator(S[beta][k])) for k in keys(S[beta])]
end