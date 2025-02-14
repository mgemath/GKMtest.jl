# example for P^3
P3 = GKMproj_space(3);
R = equivariant_cohomology_ring(P3);
con = build_GKM_connection(P3);

# build decorated tree, which is 1-2-3, mapping 1 -> 1, 2-> 4, 3->3, and marked points (1, 2) at 1, and (3) at 2.
tree = Graph{Undirected}(2);
add_edge!(tree, 1, 2);
#add_edge!(tree, 2, 3);
#add_edge!(tree, 1, 3);

#vDict = [1, 4, 3];
vDict = [1, 4]

# edgeMult = Dict(Edge(1, 2) => 1, Edge(2, 3) => 1);
#edgeMult = Dict(Edge(1, 2) => 1, Edge(1, 3) => 1);
edgeMult = Dict(Edge(1, 2) => 1)

marks = [1, 2]

dt = decoratedTree(P3, tree, vDict, edgeMult, marks, R);

C = R.coeffRing;
(t1, t2, t3, t4) = gens(C);
P = ev(1, pointClass(1, R)) * ev(2, pointClass(4, R))

println("Contribution of Edge(1,4) with d=1 in P3:")
println(GWTreeContribution(dt, R, con, P))


println("Contribution of Edge(1,4) with d=2 in P3:")
edgeMult = Dict(Edge(1, 2) => 1)
dt = decoratedTree(P3, tree, vDict, edgeMult, marks, R);
println(GWTreeContribution(dt, R, con, P))


#Warning: the example below does not work because vertices(G) throws an error for graphs G without edges.
#Conclusion: We always need special treetment for beta = zero.
println("Contribution of single point as decorated tree:")
marks = [1, 1]
tree = Graph{Undirected}(1)
vDict = [1]
dt = decoratedTree(P3, tree, vDict, edgeMult, marks, R);
println(GWTreeContribution(dt, R, con, P))