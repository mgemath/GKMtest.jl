"""
Build the second homology of the GKM graph.
Note that the returned ans.H2 is by definition the quotient of the Z-module generated by the edges
modulo the relations in second integral homology. Thus, ans.H2 is the submodule of second integral homology
generated by classes of T-invariant curves.
"""
function GKM_second_homology(G::AbstractGKM_graph)::GKM_H2

  nEdges = n_edges(G.g)
  edgeList = Vector{Edge}()
  sizehint!(edgeList, nEdges)
  for e in edges(G.g)
    push!(edgeList, e)
  end

  edgeToGenIndex = Dict{Edge, Int64}()
  for i in 1:length(edgeList)
    e = edgeList[i]
    edgeToGenIndex[e] = i
    edgeToGenIndex[reverse(e)] = i
  end

  r = rank_torus(G)
  M = free_module(ZZ, nEdges)

  cycles = _calculate_graph_cycles(G, edgeList, M)
  relations = Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  sizehint!(relations, r * length(cycles))
  
  for c in cycles
    for i in 1:r
      rel = zero(M)
      for (e,mult) in c
        edgeIndex = indexin([e], edgeList)[1]
        rel += mult * G.w[e][i] * gens(M)[edgeIndex]
      end
      push!(relations, rel)
    end
  end

  R, _ = sub(M, relations)
  H2withTorsion, q1 = quo(M, R)
  H2, q2 = _remove_torsion(H2withTorsion)

  edgeLattice = M
  quotientMap = compose(q1, q2)

  # generate dual cone of effective cone in H2.
  rkH2 = length(gens(H2))
  ImodElts = [ quotientMap(gens(edgeLattice)[i]) for i in 1:nEdges]
  I = [[v[j] for j in 1:rkH2] for v in ImodElts]
  C = cone_from_inequalities(-I) # this is the dual cone of the T-invariant curve classes
  s = sum(rays(C))
  
  # normalize s so that the minimum of edgeCurveClasses evaluated on s is 1.
  minEval = QQ(-1)
  for e in edges(G.g)
    eClass = quotientMap(gens(edgeLattice)[edgeToGenIndex[e]])
    se = sum([s[i] * eClass[i] for i in 1:rkH2])

    @req se > 0 "Edge curve class evaluates negatively on positive cone element s!"

    if se < minEval || minEval < 0
      minEval = se
    end
  end
  dualConeRaySum = (1//minEval) * s

  return GKM_H2(G, edgeLattice, H2, edgeToGenIndex, quotientMap, dualConeRaySum)
end

"""
Return a basis of the first homology of graph underlying the GKM graph.
"""
function _calculate_graph_cycles(G::AbstractGKM_graph, edgeList::Vector{Edge}, M::AbstractAlgebra.Generic.FreeModule{ZZRingElem})::Vector{Vector{Tuple{Edge, ZZRingElem}}}

  nVertices = n_vertices(G.g)
  nEdges = n_edges(G.g)

  C0 = free_module(ZZ, nVertices)
  C1 = M

  dVals = Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  sizehint!(dVals, nEdges)

  for e in edgeList
    de = gens(C0)[src(e)] - gens(C0)[dst(e)]
    push!(dVals, de)
  end
  
  d = ModuleHomomorphism(C1, C0, dVals)
  K, k = kernel(d)

  cycles = Vector{Vector{Tuple{Edge, ZZRingElem}}}()
  for cyc in gens(K)
    cycEdges = Vector{Tuple{Edge, ZZRingElem}}()
    cycIm = k(cyc)
    for i in 1:nEdges
      push!(cycEdges, (edgeList[i], cycIm[i]))
    end
    push!(cycles, cycEdges)
  end

  return cycles
end

function _remove_torsion(M::AbstractAlgebra.FPModule{ZZRingElem})

  normalForm, s = snf(M)
  invariantFactors = normalForm.invariant_factors
  torsionGenerators = Vector{}()

  for i in 1:length(invariantFactors)
    if invariantFactors[i] != 0
      push!(torsionGenerators, s(gens(normalForm)[i]))
    end
  end

  torsionSub, _ = sub(M, torsionGenerators)

  return quo(M, torsionSub)
end

function Base.show(io::IO, h2::GKM_H2)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM curve classes in H_2")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM curve classes in H_2 for $(h2.gkm)")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain",  h2::GKM_H2)
  print(io, "Curve classes in H_2 for $(h2.gkm): $(h2.H2)")
end

"""
  Return the second homology class represented by the given edge.
"""
function edgeCurveClass(H2::GKM_H2, e::Edge)
  i = H2.edgeToGenIndex[e]
  return H2.quotientMap(gens(H2.edgeLattice)[i])
end

"""
Return an iterator over the multiplicities for the given edges that sum to the second homology class beta.
"""
function _multiplicities(
  H2::GKM_H2,
  edges::Vector{Edge},
  beta::AbstractAlgebra.Generic.QuotientModuleElem{ZZRingElem};
  check::Bool = true
)
  if check
    @req parent(beta) == H2.H2 "Beta must be element of H2.H2"
  end

  nTreeEdges = length(edges)
  T = free_module(ZZ, nTreeEdges)
  t = ModuleHomomorphism(T, H2.edgeLattice, [gens(H2.edgeLattice)[H2.edgeToGenIndex[e]] for e in edges])
  q = compose(t, H2.quotientMap)

  e0 = nothing
  try
    e0 = preimage(q, beta)
    #println("e0: $e0")
  catch err
    if isa(err, ArgumentError)
      # In this case, beta cannot be expressed in terms of the given edges.
      #println("beta not in span of edge classes.")
      return Vector{}()
    else
      rethrow(err)
    end
  end

  K, k = kernel(q)
  rk = rank(K)
  mk = transpose(matrix(k)) # columns are images of generators
  #println("mk: $mk")

  P = polyhedron(-mk, [e0[i] for i in 1:nTreeEdges])
  lpts = interior_lattice_points(P)

  return (e0 + k(K([v[i] for i in 1:rk])) for v in lpts)
end

function all_classes(G::AbstractGKM_graph)
  H2 = GKM_second_homology(G);

  return [edgeCurveClass(H2, e) for e in edges(G.g)]
  # d = Dict()
  # for e in edges(G.g)
  #   d[]
  # end
end

"""
Return an upper bound for the number of edges that can be used to represent beta.
This is calculated by taking inner product with the sum of the rays of the dual cone of all edge curve classes.
If beta is not representable as positive sum of edge curve classes then this might be negative.
"""
function _max_n_edges(H2::GKM_H2, beta::AbstractAlgebra.Generic.QuotientModuleElem{ZZRingElem})
  rkH2 = length(gens(H2.H2))
  s = sum([H2.dualConeRaySum[i] * beta[i] for i in 1:rkH2])
  return Int64(ZZ(floor(s)))
end