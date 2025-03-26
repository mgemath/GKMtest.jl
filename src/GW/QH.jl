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

  # Generating P_input takes a little if we have many vertices. So we calculate it here and reuse it.
  println("Building classes to integrate over the moduli space...")
  nv = n_vertices(G.g)
  evClasses = [ev(i, pointClass(j, G)) for i in 1:3, j in 1:nv]
  P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
  println("Starting integration:")

  for c in 0:maxChernNumber
    for beta in _effectiveClassesWithChernNumber(GKM_second_homology(G), c)
      println("Calculating structure constants for $beta, Chern number $c:")
      QH_structure_constants(G, beta; refresh, P_input)
    end
  end

  G.know_all_QH_structure_consts = true

  return G.QH_structure_consts
end

function QH_structure_constants(G::AbstractGKM_graph, beta::CurveClass_type; refresh::Bool=false, P_input=nothing)

  if refresh || !haskey(G.QH_structure_consts, beta)

    R = G.equivariantCohomology.coeffRingLocalized
    nv = n_vertices(G.g)
    res = zeros(R, nv, nv, nv)
  
    if beta == zero(parent(beta))
      for i in 1:nv
        res[i, i, i] = (eulerClass(i, G)//1)^2
      end
      G.QH_structure_consts[beta] = res
      return res
    end
    if !isEffectiveCurveClass(G, beta)
      G.QH_structure_consts[beta] = res
      return res
    end
  
    # If we calculate for multiple beta, it makes sense to reuse P_input, hence the optional argument.
    # Generating P_input actually takes a surprisingly long amount of time for growing number of vertices.
    # For example for n=16 vertices it takes around 10 seconds.
    if isnothing(P_input)
      println("Building classes to integrate over the moduli space...")
      evClasses = [ev(i, pointClass(j, G)) for i in 1:3, j in 1:nv]
      P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
      println("Starting integration:")
    end
    res = integrateGKM(G, beta, 3, P_input)
    G.QH_structure_consts[beta] = res
    return res
  end
  return G.QH_structure_consts[beta]
end

function quantumProduct(
  G::AbstractGKM_graph,
  beta::CurveClass_type,
  class1, #class1 and class2 should be free module elements over the coefficient ring (or its frac field)
  class2;
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

function QH_Structure_constants_in_basis(G::AbstractGKM_graph, b::Matrix)
  s = size(b)
  @req s[1] == s[2] "Base matrix must be square"
  @req s[1] == n_vertices(G.g) "dimension of basis elements must be number of vertices of GKM graph"

  FF = G.equivariantCohomology.coeffRingLocalized
  nv = n_vertices(G.g)
  bMatrix = zero_matrix(FF, nv, nv)
  for k in keys(b)
    bMatrix[k] = FF(b[k])
  end
  @req is_invertible(bMatrix) "Base matrix must be invertible over fraction field"
  bMatrixInv = inv(transpose(bMatrix))

  res = Dict{CurveClass_type, Array{Any, 3}}()
  g = gens(G.equivariantCohomology.cohomRingLocalized)
  baseClasses = [sum([bMatrix[i,j] * g[j] for j in 1:nv]) for i in 1:nv]
  for beta in keys(G.QH_structure_consts)
    resBeta = zeros(G.equivariantCohomology.coeffRingLocalized, (nv, nv, nv)...)
    for i in 1:nv, j in 1:nv
      prodIJ = quantumProduct(G, beta, baseClasses[i], baseClasses[j])
      prodIJVect = [prodIJ[k] for k in 1:nv]
      inBasis = bMatrixInv * prodIJVect
      resBeta[i,j,:] = inBasis
    end
    res[beta] = resBeta
  end
  return res
end

function quantum_product_at_q1(G::AbstractGKM_graph, class)

  SC = QH_structure_constants(G)
  nv = n_vertices(G.g)

  M = zero_matrix(G.equivariantCohomology.coeffRingLocalized, nv, nv)
  g = gens(G.equivariantCohomology.cohomRingLocalized)

  for beta in keys(SC)
    for i in 1:nv
      cg = quantumProduct(G, beta, class, g[i])
      M[i,:] += [cg[j] for j in 1:nv]

    end
  end

  return M
end

function c1_at_q1(G::AbstractGKM_graph)
  return quantum_product_at_q1(G, firstChernClass(G))
end

"""
Return the eigenvalues of quantum multiplication by the first Chern class of the tangent bundle at q=1, t=0,
where t are the equivariant parameters.
"""
function conjecture_O_eigenvalues(G::AbstractGKM_graph; printData::Bool=true)
  c1Mat = c1_at_q1(G)
  chi = charpoly(c1Mat)
  chi0 = polynomial(QQ, [0])
  z = repeat([0], rank_torus(G))
  for i in 0:(length(chi)-1)
    set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  end
  printData && println("Characteristic poly of c1(TX)* at q=1, t=0:\n$chi0")
  return roots(QQBar, chi0)
end

function QH_is_associative(G::AbstractGKM_graph; printDiagnostics::Bool = true)::Bool
  nv = n_vertices(G.g)
  g = [QH_class(G, x) for x in gens(G.equivariantCohomology.cohomRingLocalized)]
  for i in 1:nv, j in 1:nv, k in 1:nv
    if (g[i]*g[j])*g[k] != g[i]*(g[j]*g[k])
      printDiagnostics && println("Not associative for: ($i, $j, $k)")
      return false
    end
  end
  return true
end

function QH_is_commutative(G::AbstractGKM_graph)::Bool
  nv = n_vertices(G.g)
  g = [QH_class(G, x) for x in gens(G.equivariantCohomology.cohomRingLocalized)]
  for i in 1:nv, j in 1:nv
    if g[i] * g[j] != g[j] * g[i]
      return false
    end
  end
  return true
end