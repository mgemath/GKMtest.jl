# This is dedicated to computations using P1

P1 = GKMproj_space(1);
R = equivariant_cohomology_ring(P1);
C= R.coeffRing;
(t1, t2) = gens(C);

c1 = pointClass(1, R); #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(2, R);

P = class_one();
@test integrateGKM(P1, 0, P, 1) == 1 #expected 1

P = ev(1, c1);
@test integrateGKM(P1, 1, P, 1) == 1

P = ev(1, c1)*ev(2, c1)*ev(3, c2)
@test integrateGKM(P1, 3, P, 1) == 1

@test integrateGKM(P1, 0, class_one(), 3) == 0
@test integrateGKM(P1, 2, class_one(), 3) == 0

P = Psi(1)
@test integrateGKM(P1, 1, P, 1) == -2

P = Psi(1, 1)
@test integrateGKM(P1, 2, P, 1) == 2
P = Psi(2, 2)
@test integrateGKM(P1, 2, P, 2) == 5//4