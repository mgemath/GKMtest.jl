function Seidel_element(G::AbstractGKM_graph,
    w::AbstractAlgebra.Generic.FreeModuleElem{R},
    cMax::Int64
  ) where R<:GKM_weight_type

  @req isStrictlyNEF(G) "G is required to be strictly NEF for the full Seidel element."
  SG = Seidel_space(G, w) # by default, we take the basepoint 1 here.

  nv = n_vertices(G.g)

  # calculate minimum Chern number of v_0 -> v_inf for v a vertex of G:
  cMin = ZZ(0)
  for v in 1:nv
    cMin = min(cMin, chernNumber(Edge(v, v+nv), SG))
  end

  #initialize classes to integrate over:
  res = Dict{CurveClass_type, Any}()
  P = [ev(1, point_class(nv+i, SG)) for i in 1:nv]
  coeffRing = G.equivariantCohomology.coeffRingLocalized
  g = gens(coeffRing)
  SGCC = GKM_second_homology(SG)
  GH2proj = SGCC.verticalProjection

  for c in cMin:cMax
    for b in _effectiveSectionClassesWithChernNumber(SG, c)
      println("Chern number $c, curve class $b:")
      GW = integrateGKM(SG, b, 1, P; check_degrees=true) #TODO: remove this and have global flag for checking
      @req !haskey(res, b) "Cannot have two distinct section classes restricting to the same vertical curve class!"
      res[GH2proj(b)] = [evaluate(GW[i], vcat([g[j] for j in 1:nv], [zero(coeffRing)])) for i in 1:nv]
    end
  end

  return QH_class(G, res)
end