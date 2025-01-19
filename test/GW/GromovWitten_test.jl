# example for P^3
P3 = GKMproj_space(3);
R = equivariant_cohomology_ring(P3);
con = build_GKM_connection(P3);

# build decorated tree, which is 1-2-3, mapping 1 -> 1, 2-> 4, 3->3, and marked points (1, 2) at 1, and (3) at 2.
tree = Graph{Undirected}(3);
add_edge!(tree, 1, 2);
# add_edge!(tree, 2, 3);
add_edge!(tree, 1, 3);

vDict = Dict(1=>1, 2=>4, 3=>3);

# edgeMult = Dict(Edge(1, 2) => 1, Edge(2, 3) => 1);
edgeMult = Dict(Edge(1, 2) => 1, Edge(1, 3) => 1);

marks = [1, 1, 2];

dt = decoratedTree(P3, tree, vDict, edgeMult, marks);

C = R.coeffRing;
(t1, t2, t3, t4) = gens(C);
H = R.cohomRing;
c1 = pointClass(1, R); #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(1, R); #gens(H)[1] # pointClass(1, R) # H([C(0), C(0), C(0), C(1)]) # gens(H)[2]
c3 = pointClass(4, R); #gens(H)[4] # pointClass(4, R)
classes = [c1, c2, c3];

GWTreeContribution(dt, R, con, classes)

#TODO: I get something different here! I get the below, multiplied by -1*(2*t1 - t3 - t4)//(t1 - t3)
#Update: This matches the documented output of GW_test.jl.

#(t1 - t2)//(t1*t2*t3 - t1*t2*t4 - t1*t3^2 + t1*t3*t4 - t2*t3*t4 + t2*t4^2 + t3^2*t4 - t3*t4^2) with 
#-1//(t1^7*t2^2*t3^2 - 2*t1^7*t2^2*t3*t4 + t1^7*t2^2*t4^2 - t1^7*t2*t3^3 + t1^7*t2*t3^2*t4 + t1^7*t2*t3*t4^2 - t1^7*t2*t4^3 + t1^7*t3^3*t4 - 2*t1^7*t3^2*t4^2 + t1^7*t3*t4^3 - t1^6*t2^3*t3^2 + 2*t1^6*t2^3*t3*t4 - t1^6*t2^3*t4^2 - t1^6*t2^2*t3^3 - t1^6*t2^2*t3^2*t4 + 5*t1^6*t2^2*t3*t4^2 - 3*t1^6*t2^2*t4^3 + 2*t1^6*t2*t3^4 + t1^6*t2*t3^3*t4 - 4*t1^6*t2*t3^2*t4^2 - 3*t1^6*t2*t3*t4^3 + 4*t1^6*t2*t4^4 - 2*t1^6*t3^4*t4 + 6*t1^6*t3^2*t4^3 - 4*t1^6*t3*t4^4 + 2*t1^5*t2^3*t3^3 - 6*t1^5*t2^3*t3*t4^2 + 4*t1^5*t2^3*t4^3 - t1^5*t2^2*t3^4 + 4*t1^5*t2^2*t3^3*t4 - 3*t1^5*t2^2*t3^2*t4^2 - 2*t1^5*t2^2*t3*t4^3 + 2*t1^5*t2^2*t4^4 - t1^5*t2*t3^5 - 5*t1^5*t2*t3^4*t4 + 3*t1^5*t2*t3^3*t4^2 + 7*t1^5*t2*t3^2*t4^3 + 2*t1^5*t2*t3*t4^4 - 6*t1^5*t2*t4^5 + t1^5*t3^5*t4 + 6*t1^5*t3^4*t4^2 - 9*t1^5*t3^3*t4^3 - 4*t1^5*t3^2*t4^4 + 6*t1^5*t3*t4^5 - t1^4*t2^3*t3^4 - 6*t1^4*t2^3*t3^3*t4 + 9*t1^4*t2^3*t3^2*t4^2 + 4*t1^4*t2^3*t3*t4^3 - 6*t1^4*t2^3*t4^4 + t1^4*t2^2*t3^5 + 3*t1^4*t2^2*t3^4*t4 - 7*t1^4*t2^2*t3^3*t4^2 + 3*t1^4*t2^2*t3^2*t4^3 - 2*t1^4*t2^2*t3*t4^4 + 2*t1^4*t2^2*t4^5 + 3*t1^4*t2*t3^5*t4 + 2*t1^4*t2*t3^4*t4^2 - 3*t1^4*t2*t3^3*t4^3 - 8*t1^4*t2*t3^2*t4^4 + 2*t1^4*t2*t3*t4^5 + 4*t1^4*t2*t4^6 - 4*t1^4*t3^5*t4^2 - 4*t1^4*t3^4*t4^3 + 16*t1^4*t3^3*t4^4 - 4*t1^4*t3^2*t4^5 - 4*t1^4*t3*t4^6 + 4*t1^3*t2^3*t3^4*t4 + 4*t1^3*t2^3*t3^3*t4^2 - 16*t1^3*t2^3*t3^2*t4^3 + 4*t1^3*t2^3*t3*t4^4 + 4*t1^3*t2^3*t4^5 - 4*t1^3*t2^2*t3^5*t4 - 2*t1^3*t2^2*t3^4*t4^2 + 8*t1^3*t2^2*t3^3*t4^3 + 3*t1^3*t2^2*t3^2*t4^4 - 2*t1^3*t2^2*t3*t4^5 - 3*t1^3*t2^2*t4^6 - 2*t1^3*t2*t3^5*t4^2 + 2*t1^3*t2*t3^4*t4^3 - 3*t1^3*t2*t3^3*t4^4 + 7*t1^3*t2*t3^2*t4^5 - 3*t1^3*t2*t3*t4^6 - t1^3*t2*t4^7 + 6*t1^3*t3^5*t4^3 - 4*t1^3*t3^4*t4^4 - 9*t1^3*t3^3*t4^5 + 6*t1^3*t3^2*t4^6 + t1^3*t3*t4^7 - 6*t1^2*t2^3*t3^4*t4^2 + 4*t1^2*t2^3*t3^3*t4^3 + 9*t1^2*t2^3*t3^2*t4^4 - 6*t1^2*t2^3*t3*t4^5 - t1^2*t2^3*t4^6 + 6*t1^2*t2^2*t3^5*t4^2 - 2*t1^2*t2^2*t3^4*t4^3 - 7*t1^2*t2^2*t3^3*t4^4 - 3*t1^2*t2^2*t3^2*t4^5 + 5*t1^2*t2^2*t3*t4^6 + t1^2*t2^2*t4^7 - 2*t1^2*t2*t3^5*t4^3 + 2*t1^2*t2*t3^4*t4^4 + 3*t1^2*t2*t3^3*t4^5 - 4*t1^2*t2*t3^2*t4^6 + t1^2*t2*t3*t4^7 - 4*t1^2*t3^5*t4^4 + 6*t1^2*t3^4*t4^5 - 2*t1^2*t3^2*t4^7 + 4*t1*t2^3*t3^4*t4^3 - 6*t1*t2^3*t3^3*t4^4 + 2*t1*t2^3*t3*t4^6 - 4*t1*t2^2*t3^5*t4^3 + 3*t1*t2^2*t3^4*t4^4 + 4*t1*t2^2*t3^3*t4^5 - t1*t2^2*t3^2*t4^6 - 2*t1*t2^2*t3*t4^7 + 3*t1*t2*t3^5*t4^4 - 5*t1*t2*t3^4*t4^5 + t1*t2*t3^3*t4^6 + t1*t2*t3^2*t4^7 + t1*t3^5*t4^5 - 2*t1*t3^4*t4^6 + t1*t3^3*t4^7 - t2^3*t3^4*t4^4 + 2*t2^3*t3^3*t4^5 - t2^3*t3^2*t4^6 + t2^2*t3^5*t4^4 - t2^2*t3^4*t4^5 - t2^2*t3^3*t4^6 + t2^2*t3^2*t4^7 - t2*t3^5*t4^5 + 2*t2*t3^4*t4^6 - t2*t3^3*t4^7) withoutc