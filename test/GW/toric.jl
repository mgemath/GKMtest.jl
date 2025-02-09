P3 = projective_space(NormalToricVariety,3);#thespaceP^3.
X = gkm_graph(domain(blow_up(P3,ideal(gens(cox_ring(P3))[1:3]))));
H2 = GKM_second_homology(X);
R = equivariant_cohomology_ring(X);

for e in edges(X.g)
    println("$(X.labels[src(e)]) -> $(X.labels[dst(e)]) => ", edgeCurveClass(H2, e))
end

# (-1, 1)
edgeCurveClass(H2, edgeFromLabels(X, "1", "2"))
edgeCurveClass(H2, edgeFromLabels(X, "1", "3"))
edgeCurveClass(H2, edgeFromLabels(X, "2", "3"))

#(1, 0)
edgeCurveClass(H2, edgeFromLabels(X, "1", "4"))


H=edgeCurveClass(H2, edgeFromLabels(X, "5", "4")); #(0,1)
# E=edgeCurveClass(H2, edgeFromLabels(X, "1", "4")) #(1,0)
E=edgeCurveClass(H2, edgeFromLabels(X, "1", "2")) #(-1,1)
d=2;#thiscanbeanynonnegativeinteger
e=-1;#thiscanbeanynonpositiveinteger
beta=d*H+e*E;
# beta = 3*E


P= prod(i -> ev(i, pointClass(1, R)), 1:(2*d+e));
println(integrateGKM(X, H2, beta, 2*d+e, P, show_bar = false));

# I got error when E = (-1,1) and d=-e.
# I got error when E = (1,0) and d=-e.
# both error are solved after adding lines 195-197 in Main