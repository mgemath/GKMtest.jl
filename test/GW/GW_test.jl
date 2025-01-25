P3 = GKMproj_space(3);
R = equivariant_cohomology_ring(P3);

C = R.coeffRing;
(t1, t2, t3, t4) = gens(C);

c1 = pointClass(1, R); #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(1, R); #gens(H)[1] # pointClass(1, R) # H([C(0), C(0), C(0), C(1)]) # gens(H)[2]
c3 = pointClass(4, R); #gens(H)[4] # pointClass(4, R)
classes = [c1, c2, c3];

integrateGKM(P3, classes)

P = ev(1, c1)*ev(2, c3)#*ev(3, c3);
integrateGKM(P3, 3, P)
integrateGKM(P3, 3, P, 1)

P2 = GKMproj_space(2);

R2 = equivariant_cohomology_ring(P2);
C2 = R2.coeffRing;
(t1, t2, t3) = gens(C2);
c_1 = pointClass(1, R2);
c_2 = pointClass(2, R2);
c_3 = pointClass(3, R2);

P = ev(1, c_1)*ev(2, c_1)*ev(3, c_1)*ev(4, c_1)*ev(5, c_1)
integrateGKM(P2, 5, P, 2)

P = ev(1, c_1)*ev(2, c_1)
integrateGKM(P2, 2, P, 1)

P = ev(1, c_1)*ev(2, c_2)
integrateGKM(P2, 2, P, 1)

(-48*t1^6 + 144*t1^5*t2 + 144*t1^5*t3 - 80*t1^4*t2^2 - 560*t1^4*t2*t3 - 80*t1^4*t3^2 - 80*t1^3*t2^3 + 560*t1^3*t2^2*t3 + 560*t1^3*t2*t3^2 - 80*t1^3*t3^3 + 144*t1^2*t2^4 - 336*t1^2*t2^3*t3 - 336*t1^2*t2^2*t3^2 - 336*t1^2*t2*t3^3 + 144*t1^2*t3^4 - 80*t1*t2^5 + 112*t1*t2^4*t3 + 112*t1*t2^3*t3^2 + 112*t1*t2^2*t3^3 + 112*t1*t2*t3^4 - 80*t1*t3^5 + 16*t2^6 - 16*t2^5*t3 - 16*t2^4*t3^2 - 16*t2^3*t3^3 - 16*t2^2*t3^4 - 16*t2*t3^5 + 16*t3^6)//(t1^6 - 3*t1^5*t2 - 3*t1^5*t3 + t1^4*t2^2 + 13*t1^4*t2*t3 + t1^4*t3^2 + 3*t1^3*t2^3 - 13*t1^3*t2^2*t3 - 13*t1^3*t2*t3^2 + 3*t1^3*t3^3 - 2*t1^2*t2^4 - t1^2*t2^3*t3 + 21*t1^2*t2^2*t3^2 - t1^2*t2*t3^3 - 2*t1^2*t3^4 + 4*t1*t2^4*t3 - 7*t1*t2^3*t3^2 - 7*t1*t2^2*t3^3 + 4*t1*t2*t3^4 - 2*t2^4*t3^2 + 5*t2^3*t3^3 - 2*t2^2*t3^4)