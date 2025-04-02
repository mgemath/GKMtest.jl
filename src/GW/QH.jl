@doc raw"""
    QH_structure_constants(G::AbstractGKM_graph; refresh::Bool=false)

Return the structure constants of the equivariant quantum cohomology $QH_T^*(X)$ where $X$ is the GKM variety realizing the GKM graph.

!!! note
    - This requires `is_strictly_nef(G)==true`, as this guarantees that there are at most finitely many curve classes $\beta$ with non-zero coefficients for $q^\beta$.
    - If `is_strictly_nef(G)==false`, use the method of `QH_structure_constants` below that specifies a specific $\beta$.
    - As this computation might be expensive, the result is stored in `G` for later use. If the requested structure constants have been computed before, they will not be
      computed afresh unless the optional argument `refresh` is set to `true`.

# Output format:
The output type is `Dict{CurveClass_type, Array{Any, 3}}`.
If `ans` denotes the returned object, then
`ans[beta][i, j, k]` is the the $q^\beta$-coefficient of $PD(v_i) \ast PD(v_j)$ localized at $v_k$, where $v_i, v_j, v_k$ represent the fixed points with indices $i,j,k$, respectively, 
and $PD$ represents the Poincaré dual.

# Optional arguments:
 - `refresh::Bool`: `false` by default. If `true`, then this will overwrite any previously calculated $QH_T$ structure constants of `G`.

# Example
```jldoctest QH_structure_constants_all
julia> P1 = projective_space(GKM_graph, 1);

julia> S = QH_structure_constants(P1; show_progress=false)
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Array{Any, 3}} with 2 entries:
  (0) => [t1^2 - 2*t1*t2 + t2^2 0; 0 0;;; 0 0; 0 t1^2 - 2*t1*t2 + t2^2]
  (1) => [1 1; 1 1;;; 1 1; 1 1]
```
"""
function QH_structure_constants(G::AbstractGKM_graph; refresh::Bool=false, show_progress::Bool=true)

  @req is_strictly_nef(G) "G is not strictly NEF, so need to specify beta"

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
  show_progress && println("Building classes to integrate over the moduli space...")
  nv = n_vertices(G.g)
  evClasses = [ev(i, point_class(j, G)) for i in 1:3, j in 1:nv]
  P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
  show_progress && println("Starting integration:")

  for c in 0:maxChernNumber
    for beta in _effectiveClassesWithChernNumber(GKM_second_homology(G), c)
      show_progress && println("Calculating structure constants for $beta, Chern number $c:")
      QH_structure_constants(G, beta; refresh = refresh, P_input = P_input, show_progress = show_progress)
    end
  end

  G.know_all_QH_structure_consts = true

  return G.QH_structure_consts
end

@doc raw"""
    QH_structure_constants(G::AbstractGKM_graph, beta::CurveClass_type; refresh::Bool=false, P_input=nothing, show_progress::Bool=true)

Return the $q^\beta$-coefficients of the structure constants of the equivariant quantum cohomology $QH_T^*(X)$, where $X$ is the GKM variety realizing the GKM graph.

!!! note
    - As this computation might be expensive, the result is stored in `G` for later use. If the requested structure constants have been computed before, they will not be
      computed afresh unless the optional argument `refresh` is set to `true`.

# Output format:
The output type is `Array{Any, 3}`.
If`ans` denotes the returned object, then
`ans[i, j, k]` is the the $q^\beta$-coefficient of $PD(v_i) \ast PD(v_j)$ localized at $v_k$, where $v_i, v_j, v_k$ represent the fixed points with indices $i,j,k$, respectively, 
and $PD$ represents the Poincaré dual.

# Optional arguments:
 - `refresh::Bool`: `false` by default. If `true`, then this will overwrite any previously calculated $QH_T$ structure constants of `G`.

# Example
```jldoctest QH_structure_constants_edge
julia> P1 = projective_space(GKM_graph, 1);

julia> beta = curve_class(P1, Edge(1, 2));

julia> QH_structure_constants(P1, 0*beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 t1^2 - 2*t1*t2 + t2^2  0
 0                      0

[:, :, 2] =
 0  0
 0  t1^2 - 2*t1*t2 + t2^2

julia> QH_structure_constants(P1, beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 1  1
 1  1

[:, :, 2] =
 1  1
 1  1

julia> QH_structure_constants(P1, 2*beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 0  0
 0  0

[:, :, 2] =
 0  0
 0  0

julia> QH_structure_constants(P1, -1 * beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 0  0
 0  0

[:, :, 2] =
 0  0
 0  0
```
"""
function QH_structure_constants(G::AbstractGKM_graph, beta::CurveClass_type; refresh::Bool=false, P_input=nothing, show_progress::Bool=true)

  if refresh || !haskey(G.QH_structure_consts, beta)

    R = G.equivariantCohomology.coeffRingLocalized
    nv = n_vertices(G.g)
    res = zeros(R, nv, nv, nv)
  
    if beta == zero(parent(beta))
      for i in 1:nv
        res[i, i, i] = (euler_class(i, G)//1)^2
      end
      G.QH_structure_consts[beta] = res
      return res
    end
    if !is_effective(G, beta)
      G.QH_structure_consts[beta] = res
      return res
    end
  
    # If we calculate for multiple beta, it makes sense to reuse P_input, hence the optional argument.
    # Generating P_input actually takes a surprisingly long amount of time for growing number of vertices.
    # For example for n=16 vertices it takes around 10 seconds.
    if isnothing(P_input)
      show_progress && println("Building classes to integrate over the moduli space...")
      evClasses = [ev(i, point_class(j, G)) for i in 1:3, j in 1:nv]
      P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
      show_progress && println("Starting integration:")
    end
    res = gromov_witten(G, beta, 3, P_input; show_bar=show_progress)
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
      eulerI = euler_class(i, G)
      eulerJ = euler_class(j, G)
      res += (C[i, j, k] * class1[i] * class2[j] // eulerI // eulerJ) * gens(G.equivariantCohomology.cohomRingLocalized)[k]
    end
    return res
  end

  P_input = [ev(1, class1)*ev(2, class2)*ev(3, point_class(v, G)) for v in 1:nv]
  GW_invts = gromov_witten(G, beta, 3, P_input)

  return sum([GW_invts[i] * gens(G.equivariantCohomology.cohomRingLocalized)[i] for i in 1:nv])
end


@doc raw"""
    QH_structure_constants_in_basis(G::AbstractGKM_graph, b::Matrix)

Return all structure constants of `G` that have been calculated so far with respect to the given basis.
A smart choice of basis can drastically simplify the presentation of the ring $QH_T^*(X)$.

!!! note
    This does not calculate any structure constants afresh.
    To do so, use `QH_structure_constants`.

# Output format:
The same as that of `QH_structure_constants`, i.e. of type `Dict{CurveClass_type, Array{Any, 3}}`.

# Arguments
 - `G::AbstractGKM_graph`: The GKM graph whose quantum cohomology is of interest.
 - `b::Matrix`: A matrix whose rows are the desired $H_T^*(\text{pt};\mathbb{Q})$-linear basis of $H_T^*(X;\mathbb{Q})$.
    The element `b[i,j]` is the localization to the `j`-th fixed point of the `i`-th basis element.

# Examples
The following example shows that $QH_T(X;\mathbb{Q}) \cong\mathbb{Q}[t_1, t_2, e]/(e^2 - (t_1-t_2)e - q)$ where $e=PD([1:0])$ and $q$
corresponds to the curve class $[\mathbb{P}^1]\in H_2(\mathbb{P}^1;\mathbb{Z})$.
```jldoctest
julia> P1 = projective_space(GKM_graph, 1);

julia> QH_structure_constants(P1; show_progress=false);

julia> P1 = projective_space(GKM_graph, 1);

julia> QH_structure_constants(P1; show_progress=false)
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Array{Any, 3}} with 2 entries:
  (0) => [t1^2 - 2*t1*t2 + t2^2 0; 0 0;;; 0 0; 0 t1^2 - 2*t1*t2 + t2^2]
  (1) => [1 1; 1 1;;; 1 1; 1 1]

julia> t1, t2 = gens(P1.equivariantCohomology.coeffRing);

julia> base = [1 1; t1-t2 0 ];

julia> QH_structure_constants_in_basis(P1, base)
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Array{Any, 3}} with 2 entries:
  (0) => [1 0; 0 0;;; 0 1; 1 t1 - t2]
  (1) => [0 0; 0 1;;; 0 0; 0 0]
```

**TODO** Add GP2 example as well!
"""

function QH_structure_constants_in_basis(G::AbstractGKM_graph, b::Matrix)
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
  return quantum_product_at_q1(G, first_chern_class(G))
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

function QH_is_polynomial(G::AbstractGKM_graph)::Bool
  S = G.QH_structure_consts
  for b in keys(S)
    for k in keys(S[b])
      s = S[b][k]
      if !_is_polynomial(s)
        println("QH is not polynomial for curve class $b at index $k.")
        return false
      end
    end
  end
  return true
end

function QH_is_homogeneous(G::AbstractGKM_graph)::Bool
  S = G.QH_structure_consts
  for b in keys(S)
    for k in keys(S[b])
      s = S[b][k]
      if !_is_homogeneous(s)
        println("QH is not homogeneous for curve class $b at index $k.")
        return false
      end
    end
  end
  return true
end

# Return all QH structure constants that are nonzero and have been calculated.
# This does not calculate potentially non-zero constants that have not yet been calculated.
function QH_supporting_curve_classes(G::AbstractGKM_graph)
  S = G.QH_structure_consts
  res = CurveClass_type[]
  for b in keys(S)
    if !all(s -> iszero(s), S[b])
      push!(res, b)
    end
  end
  return res
end