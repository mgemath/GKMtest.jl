P3 = GKMproj_space(3);
R = equivariant_cohomology_ring(P3);

(a,b,c,d)=gens(R.cohomRing)

function brute_force(G::AbstractGKM_graph, beta, tree, col; MAX = 10)

  ans = NTuple{n_edges(tree),Int64}[]
  for mul in Interators.product([1:MAX for _ in 1:n_edges(tree)]...)
    
  end
end