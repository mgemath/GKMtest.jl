G = flag_gkm_graph([1, 3]);
con = build_GKM_connection(G);


S = GKMsubgraph_from_vertices(G, [1]);
(blowupSub, blowupCon) = blowupGKM(S, con);

X = blowupSub.super;
H2 = GKM_second_homology(X);
R = equivariant_cohomology_ring(X);
H = edgeCurveClass(H2, edgeFromLabels(X, "2", "3"));
E = edgeCurveClass(H2, edgeFromLabels(X, "[1>3]", "[1>2]"));

d=2;#thiscanbeanynonnegativeinteger
e=-1;#thiscanbeanynonpositiveinteger
beta=d*H+e*E;
# beta = 3*E


P= prod(i -> ev(i, pointClass(1, R)), 1:(2*d+e));
println(integrateGKM(X, H2, beta, 2*d+e, P, show_bar = false));
