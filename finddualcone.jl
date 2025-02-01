using Oscar, ToricAtiyahBott

P3 = projective_space(NormalToricVariety, 3);# the space P^3.
Y = domain(blow_up(P3, [1,1,1]; coordinate_name="Ex1"));
X = domain(blow_up(Y, [-1,0,0]; coordinate_name="Ex2")); 

mgP3 = moment_graph(P3);
mgX = moment_graph(X);

RP3 = cohomology_ring(P3)
RX = cohomology_ring(X)

betti_number(X,2)

gens(RX)

gens(picard_group(X))
integrate(mgX[(1,2)]*cohomology_class(X, gens(cohomology_ring(X))[1]))
rational_equivalence_class(gens(RX)[1])
[-integrate(mgX[(1,2)]*cohomology_class(X, gens(RX)[i])) for i in 1:length(gens(RX))]
i=[[Int64(-integrate(mgP3[e]*cohomology_class(P3, gens(RP3)[i]))) for i in 1:length(gens(RP3))] for e in keys(mgP3)]

# remove linearly dependent divisor
# typeof(gens(RP3)[1])
res = MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}[]
for i in 1:length(gens(RP3))
end

RX
base_ring(RX)
# projectio = hom(base_ring(RX), RX, [gens(RX)[i] for i in 1:length(gens())])
projectio = hom(base_ring(RX), RX, gens(RX))
k=kernel(projectio)
gens(k)
[degree(t) == degree(gens(k)[1]) for t in gens(k)]
[t for t in gens(k) if degree(t) == degree(gens(k)[1])]
dim(base_ring(RX))

dim(base_ring(RX))
ideal([t for t in gens(k) if degree(t) == degree(gens(k)[1])])
Q, f =quo(base_ring(RX), ideal([t for t in gens(k) if degree(t) == degree(gens(k)[1])]))
Q
dim(Q)
gens(Q)
minimal_generating_set(ideal([t for t in gens(k) if degree(t) == degree(gens(k)[1])]))
s = sub(Q, gens(Q))

I = Matrix{Int64}(undef,div(length(keys(mgP3)), 2), length(gens(RP3)))
index = 1
for e in keys(mgP3)
    issorted(e) || continue
    I[index, :] = [Int64(-integrate(mgP3[e]*cohomology_class(P3, gens(RP3)[i]))) - rand(Int) for i in 1:length(gens(RP3))]
    index += 1
end

I
C=cone_from_inequalities(I)
rays(C)




[Int64(-integrate(mgP3[1,2]*cohomology_class(P3, gens(RP3)[i]))) for i in 1:length(gens(RP3))]
I[1,:] = [Int64(-integrate(mgP3[1,2]*cohomology_class(P3, gens(RP3)[i]))) for i in 1:length(gens(RP3))]
foreach(i -> I[], div(keys(mgP3), 2))
I
i
push!(I, i[1])
cone_from_inequalities

typeof([1 0; 1 1])

matrix([[1,0], [1,1]]) == [1 0; 1 1]
typeof(matrix([[1,0], [1,1]]))

