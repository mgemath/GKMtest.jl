"""
Build the second homology of the GKM graph.
Note that the returned ans.H2 is by definition the quotient of the Z-module generated by the edges
modulo the relations in second integral homology. Thus, ans.H2 is the submodule of second integral homology
generated by classes of T-invariant curves.

Note: If G.curveClasses is already set, it just returns that. Otherwise it calculates it afresh and stores it
there.
"""
function GKM_second_homology(G::AbstractGKM_graph)::GKM_H2
  if isnothing(G.curveClasses)
    G.curveClasses = _GKM_second_homology(G)
  end
  return G.curveClasses
end


function _GKM_second_homology(G::AbstractGKM_graph)::GKM_H2

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
  cwd = _common_weight_denominator(G)

  for c in cycles
    for i in 1:r
      rel = zero(M)
      for (e,mult) in c
        edgeIndex = indexin([e], edgeList)[1]
        rel += mult * ZZ(cwd * G.w[e][i]) * gens(M)[edgeIndex]
      end
      push!(relations, rel)
    end
  end

  R, _ = sub(M, relations)
  H2withTorsion, q1 = quo(M, R)
  H2, q2 = _remove_torsion(H2withTorsion)

  edgeLattice = M
  quotientMap = compose(q1, q2)

  dualConeRaySum, C, H2ToCN = _finish_GKM_H2(edgeLattice, H2, quotientMap, G, edgeToGenIndex)

  return GKM_H2(G, edgeLattice, H2, edgeToGenIndex, quotientMap, dualConeRaySum, C, H2ToCN, nothing)
end

function _finish_GKM_H2(edgeLattice, H2, quotientMap, G, edgeToGenIndex)

  nEdges = n_edges(G.g)

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

  # generate Chern numbers as ZZ-module homorphism from H2 to ZZ
  ZZasModule = free_module(ZZ, 1)
  edgeCNDict = Dict{Int64, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  edgeCNVect = Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for e in edges(G.g)
    edgeCNDict[edgeToGenIndex[e]] = chernNumber(e, G) * gens(ZZasModule)[1]
  end
  for i in 1:nEdges
    push!(edgeCNVect, edgeCNDict[i])
  end
  edgeLatticeToCN = ModuleHomomorphism(edgeLattice, ZZasModule, edgeCNVect)

  H2ToCN = Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for g in gens(H2)
    p = preimage(quotientMap, g)
    push!(H2ToCN, edgeLatticeToCN(p))
  end
  H2ToCN = ModuleHomomorphism(H2, ZZasModule, H2ToCN)

  return (dualConeRaySum, C, H2ToCN)
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
function edgeCurveClass(G::AbstractGKM_graph, src::String, dst::String)
  @req src in G.labels "Source vertex not found for label $src"
  @req dst in G.labels "Destination vertex not found for label $dst"
  sd = indexin([src, dst], G.labels)
  return edgeCurveClass(G, Edge(sd[1], sd[2]))
end

"""
  Return the second homology class represented by the given edge.
"""
function edgeCurveClass(G::AbstractGKM_graph, e::Edge)
  return edgeCurveClass(GKM_second_homology(G), e)
end

"""
  Return the second homology class represented by the given edge.
"""
function edgeCurveClass(H2::GKM_H2, e::Edge)
  i = H2.edgeToGenIndex[e]
  return H2.quotientMap(gens(H2.edgeLattice)[i])
end

"""
Return whether beta is in the GW effective cone, i.e. whether it is a non-negative
linear combination of edge curve classes. Note that we consider the effective cone to be 
closed and hence zero is also considered effective.
"""
function isEffectiveCurveClass(
  G::AbstractGKM_graph,
  beta::CurveClass_type
)::Bool
  # Note: It wouldn't make sense to calculate G.curveClasses afresh here, since the user
  # must have gotten beta from somewhere, and it must be compatible with G.curveClasses.
  @req !isnothing(G.curveClasses) "G.curveClasses has not been calculated yet. Where did you get beta from?"
  return isEffectiveCurveClass(G.curveClasses, beta)
end

"""
Return whether beta is in the GW effective cone, i.e. whether it is a non-negative
linear combination of edge curve classes. Note that we consider the effective cone to be 
closed and hence zero is also considered effective.
"""
function isEffectiveCurveClass(
  H2::GKM_H2,
  beta::CurveClass_type
)::Bool
  rkH2 = length(gens(H2.H2))
  return all(r -> sum([beta[i] * r[i] for i in 1:rkH2 ]) >= 0, rays(H2.dualCone))
end

"""
Return the chern number of the curve class represented by the given edge.
This is the pairing of the curve class with the first Chern class of the tangent bundle.
"""
function chernNumber(e::Edge, G::AbstractGKM_graph)::ZZRingElem

  R = G.equivariantCohomology
  cn = integrateOverEdge(firstChernClass(R), R, e)
  (flag, quotient) = divides(numerator(cn), denominator(cn))

  @req flag "1st Chern class not a GKM class!"
  @req is_constant(quotient) "1st Chern number has too high degree"

  if length(coefficients(quotient)) > 0
    return ZZ(coeff(quotient, 1))
  else
    return ZZ(0)
  end
end

"""
Return the chern number of the second homology class.
This is the pairing of the second homology class with the first Chern class of the tangent bundle.
"""
function chernNumber(G::AbstractGKM_graph, beta::CurveClass_type; check::Bool=true)::ZZRingElem
  # Note: It wouldn't make sense to calculate G.curveClasses afresh here, since the user
  # must have gotten beta from somewhere, and it must be compatible with G.curveClasses.
  @req !isnothing(G.curveClasses) "G.curveClasses has not been calculated yet. Where did you get beta from?"
  return chernNumber(G.curveClasses, beta; check)
end

"""
Return the chern number of the second homology class.
This is the pairing of the second homology class with the first Chern class of the tangent bundle.
"""
function chernNumber(H2::GKM_H2, beta::CurveClass_type; check::Bool=true)::ZZRingElem
  if check
    @req parent(beta) == H2.H2 "Beta must be element of H2.H2"
  end
  return H2.chernNumber(beta)[1]
end

"""
Return true if and only if the Chern numbers of all curve classes corresponding to
edges of the GKM graph are strictly positive.
"""
function isStrictlyNEF(G::AbstractGKM_graph)::Bool
  for e in edges(G.g)
    if chernNumber(e, G) <= 0
      return false
    end
  end
  return true
end

"""
Return all effective classes in H2.H2 with prescribed Chern number (i.e. class evaluated on c1(TX))
Note that this considers the zero class as effective.

Warning: This only works if isStrictlyNEF(H2.gkm, R) holds, as otherwise there might be infinitely
many effective curve classes of given Chern number.
"""
function _effectiveClassesWithChernNumber(
  H2::GKM_H2,
  chernNumber::Int64;
  check::Bool=true
)

  if check
    @req isStrictlyNEF(H2.gkm) "Space is not strictly nef, so there are infinitely many effective curve classes with bounded chern number."
  end

  q = H2.chernNumber

  ZZasModule = codomain(H2.chernNumber)
  ZZgen = gens(ZZasModule)[1]
  e0 = nothing
  try
    e0 = preimage(H2.chernNumber, chernNumber * ZZgen)
  catch err
    if isa(err, ArgumentError)
      # In this case, no integer combination of edge curve classes has this chern number.
      return Vector{}()
    else
      rethrow(err)
    end
  end

  K, k = kernel(q)
  rk = rank(K)
  mk = transpose(matrix(k)) # columns are images of generators

  dualRays = rays(H2.dualCone)
  nDualRays = length(dualRays)
  rkH2 = length(gens(H2.H2))
  rayMatrix = QQMatrix(nDualRays, rkH2)
  for i in 1:nDualRays
    for j in 1:rkH2
      rayMatrix[i, j] = dualRays[i][j]
    end
  end
  
  Re0 = rayMatrix*[e0[i] for i in 1:rkH2]
  P = polyhedron(-rayMatrix*mk, Re0)

  ptsIterator = (e0 + k(K([v[i] for i in 1:rk])) for v in lattice_points(P))
  return ptsIterator
end

"""
Return an iterator over the multiplicities for the given edges that sum to the second homology class beta.
Each entry in each returned multiplicity is positive.
Each returned multiplicity has type Vector{Int64}.
"""
function _multiplicities(
  H2::GKM_H2,
  edges::Vector{Edge},
  beta::CurveClass_type;
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

  P = polyhedron(-mk, [e0[i] for i in 1:nTreeEdges])
  lpts = interior_lattice_points(P)

  intPtsIterator = (e0 + k(K([v[i] for i in 1:rk])) for v in lpts)

  # Make sure all multiplicity-sets have positive entries.
  # The reason for potentially zero entries here is that Oscar's "interior_lattice_points()" works intrinsically 
  # with respect to the polytope and not with respect to the defining inequalities.
  # For example, the polytope in R^2 defined by 0 <= x <= 0 and -2 <= y <= 2 has interior lattice points
  # (0,-1), (0, 0), (0, 1), even though the first defining inequality is not strict at these points.
  noZerosIterator = Iterators.filter(m -> all(i -> m[i] > 0, 1:nTreeEdges), intPtsIterator)

  # Convert to Vect{Int64}:
  # If we get an InexactError here, it means that m[i] is too large to be converted ti Int64.
  # This is extremely unlikely in the context of this package.
  return Iterators.map(m -> [Int64(m[i]) for i in 1:nTreeEdges], noZerosIterator)
end

function all_classes(G::AbstractGKM_graph)
  H2 = GKM_second_homology(G)

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
function _max_n_edges(H2::GKM_H2, beta::CurveClass_type)
  rkH2 = length(gens(H2.H2))
  s = sum([H2.dualConeRaySum[i] * beta[i] for i in 1:rkH2])
  return Int64(ZZ(floor(s)))
end