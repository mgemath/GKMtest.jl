# This is dedicated to computations using P1

P2 = GKMproj_space(2);
R = equivariant_cohomology_ring(P2);
C= R.coeffRing;
(t1, t2, t3) = gens(C);

c1 = pointClass(1, R); #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(2, R);
c3 = pointClass(3, R);

P = ev(1, c1)*ev(2, c1);
@test integrateGKM(P2, 2, P, 1) == 1

P = ev(1, c1)*ev(2, c1)*ev(3, c1)*ev(4, c1)*ev(5, c1);
@test integrateGKM(P2, 5, P, 2) == 1

P = class_one();
@test integrateGKM(P2, 2, P, 1) == 0
@test integrateGKM(P2, 0, P, 1) == 0

P = ev(1, c1)*ev(2, c1)*ev(3, c1)*ev(4, c1)*ev(5, c1)*ev(6, c1)*ev(7, c1)*ev(8, c1);
@test integrateGKM(P2, 8, P, 3) == 12