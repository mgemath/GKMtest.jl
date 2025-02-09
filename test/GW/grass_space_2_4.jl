F = flag_gkm_graph([2,2]);

R = equivariant_cohomology_ring(F);
C= R.coeffRing;
(t1, t2, t3, t4) = gens(C);
h2 = GKM_second_homology(F);

c1 = pointClass(1, R);

# P = class_one();
# @test integrateGKM(F, 0, P, 1) == 0 #expected 0

# P = class_one();
# @test integrateGKM(F, 1, P, 1) == 0 #expected 0

P = ev(1, c1)*Psi(2);
x=integrateGKM(F, h2, 1*gens(h2.H2)[1], 1, P) #expected ?

println(x)

P = ev(1, c1)*Psi(2, 1);
x=integrateGKM(F, h2, 1*gens(h2.H2)[1], 1, P) #expected ?
