function QH_structure_constants(G::AbstractGKM_graph; refresh::Bool=false)

  @req isStrictlyNEF(G) "G is not strictly NEF, so need to specify beta"

  if !refresh && G.know_all_QH_structure_consts
    return G.QH_structure_consts
  end

  if refresh
    empty!(G.QH_structure_consts)
  end

  # Max deg of cohomology class (primitive wrt equivariant parameters): 2n
  # Max deg of product of two such: 4n
  # Degree of q^beta: 2c1(beta)
  # Thus, assuming G is strictly NEF, we only need 0 <= c1(beta) <= 2n
  n = valency(G)
  maxChernNumber = 2*n

  for c in 0:maxChernNumber
    for beta in _effectiveClassesWithChernNumber(GKM_second_homology(G), c)
      println("Calculate structure constants for $beta, Chern number $c:")
      QH_structure_constants(G, beta; refresh)
    end
  end

  G.know_all_QH_structure_consts = true

  return G.QH_structure_consts
end

function QH_structure_constants(G::AbstractGKM_graph, beta::CurveClass_type; refresh::Bool=false)

  if refresh || !haskey(G.QH_structure_consts, beta)

    R = G.equivariantCohomology.coeffRingLocalized
    nv = n_vertices(G.g)
    res = zeros(R, nv, nv, nv)
  
    if beta == zero(parent(beta))
      for i in 1:nv
        res[i, i, i] = one(R)
      end
      G.QH_structure_consts[beta] = res
      return res
    end
    if !isEffectiveCurveClass(G, beta)
      G.QH_structure_consts[beta] = res
      return res
    end
  
    P_input = [ev(1, pointClass(i, G))*ev(2, pointClass(j, G))*ev(3, pointClass(k, G)) for i in 1:nv, j in 1:nv, k in 1:nv]
    res = integrateGKM(G, beta, 3, P_input)
    G.QH_structure_consts[beta] = res
    return res
  end
  return G.QH_structure_consts[beta]
end

function quantumProduct(
  G::AbstractGKM_graph,
  beta::CurveClass_type,
  class1::FreeModElem{QQMPolyRingElem},
  class2::FreeModElem{QQMPolyRingElem};
  useStructureConstants::Bool = true
)
  if beta == 0
    return class1 * class2
  end

  nv = n_vertices(G.g)

  if useStructureConstants
    C = QH_structure_constants(G, beta)
    res = zero(G.equivariantCohomology.cohomRingLocalized)
    for i in 1:nv, j in 1:nv, k in 1:nv
      eulerI = eulerClass(i, G)
      eulerJ = eulerClass(j, G)
      res += (C[i, j, k] * class1[i] * class2[j] // eulerI // eulerJ) * gens(G.equivariantCohomology.cohomRingLocalized)[k]
    end
    return res
  end

  P_input = [ev(1, class1)*ev(2, class2)*ev(3, pointClass(v, G)) for v in 1:nv]
  GW_invts = integrateGKM(G, beta, 3, P_input)

  return sum([GW_invts[i] * gens(G.equivariantCohomology.cohomRingLocalized)[i] for i in 1:nv])
end