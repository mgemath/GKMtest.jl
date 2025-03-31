#TODO: If calculating GKM_second_homology of the Seidel space takes too long, implement it here.

@doc raw"""
    Seidel_space(G::AbstractGKM_graph, w::AbstractAlgebra.Generic.FreeModuleElem{R}) -> AbstractGKM_graph

It constructs the Seidel space of the `G` with weights `w`.
"""
function Seidel_space(
  G::AbstractGKM_graph,
  weight::AbstractAlgebra.Generic.FreeModuleElem{R};
  basePoint::Int64 = 1
) where R<:GKM_weight_type

  @req parent(weight) === G.M "The weight does not belong to the right character lattice"
  @req is_connected(G.g) "The GKM graph should be connected."

  nv = n_vertices(G.g)
  ne = n_edges(G.g)
  r = rank_torus(G)
  labels = G.labels
  
  # create labels for Seidel space (vertices 1...n over 0, n+1...2n over inf)
  Slabels = Vector{String}()
  sizehint!(Slabels, 2*nv)
  for l in labels
    push!(Slabels, "[" * l * "]_0")
  end
  for l in labels
    push!(Slabels, "[" * l * "]_inf")
  end

  # new GKM graph:
  SM = free_module(base_ring(G.M), r+1)
  Sw = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{G.weightType}}()
  SG = gkm_graph(Graph{Undirected}(2*nv), Slabels, SM, Sw; checkLabels = false)

  # weights:
  z = gens(SM)[r+1]
  iM = ModuleHomomorphism(G.M, SM, [gens(SM)[i] for i in 1:r])

  # add edges:
  for i in 1:nv
    add_edge!(SG, i, i+nv, z)
  end
  for e in edges(G.g)
    # edge over 0:
    add_edge!(SG, src(e), dst(e), iM(G.w[e]))
    # edge over inf:
    wz = sum([G.w[e][i] * weight[i] for i in 1:r])
    add_edge!(SG, nv + src(e), nv + dst(e), iM(G.w[e]) -  wz * z)
  end

  # infer GKM connection to Seidel space if possible
  con = get_connection(G)
  if !isnothing(con)
    SconDict = Dict{Tuple{Edge, Edge}, Edge}()
    # connection along new edge
    for v in 1:nv
      e = Edge(v, nv+v)
      SconDict[(e, e)] = reverse(e)
      SconDict[(reverse(e), reverse(e))] = e
      for w in all_neighbors(G.g, v)
        ei = Edge(v, w)
        epi = Edge(nv+v, nv+w)
        SconDict[(e, ei)] = epi
        SconDict[(reverse(e), epi)] = ei
      end
    end
    # connection over 0 and infty
    for v in 1:nv
      for w in all_neighbors(G.g, v)
        for u in all_neighbors(G.g, v)
          # over zero:
          e = Edge(v, w)
          ei = Edge(v, u)
          epi = con.con[(e, ei)]
          SconDict[(e, ei)] = epi
          # over infinity:
          e = Edge(v + nv, w + nv)
          ei = Edge(v + nv, u + nv)
          epi = Edge(src(epi) + nv, dst(epi) + nv)
          SconDict[(e, ei)] = epi
        end
      end
    end
    for e in edges(G.g)
      # over 0
      v = src(e)
      w = dst(e)
      ei = Edge(v, nv+v)
      epi = Edge(w, nv+w)
      SconDict[(e, ei)] = epi
      SconDict[(reverse(e), epi)] = ei
      # over infty
      e = Edge(v + nv, w + nv)
      ei = reverse(ei)
      epi = reverse(epi)
      SconDict[(e, ei)] = epi
      SconDict[(reverse(e), epi)] = ei
    end
    Scon = build_GKM_connection(SG, SconDict)
    set_connection!(SG, Scon)
  end


  # build curve classes:

  GCC = GKM_second_homology(G)
  Gq = GCC.quotientMap
  rkGH2 = length(gens(GCC.H2))
  horizontalEdgeLattice = free_module(ZZ, nv)
  # first edges over 0, then over inf, then horizontal ones.
  SedgeLattice, _, _ = direct_sum(vcat(repeat([GCC.edgeLattice], 2), [horizontalEdgeLattice]))
  extraFactor = free_module(ZZ, 1)
  SH2, _, projs = direct_sum(GCC.H2, extraFactor)

  #build map from horizontalEdgeLattice to GCC.edgeLattice and extraFactor
  horToExtra = ModuleHomomorphism(horizontalEdgeLattice, extraFactor, [gens(extraFactor)[1] for i in 1:nv])
  # calculate H2 shifts of horizontal edges
  shiftDict = Dict{Int64, CurveClass_type}()
  shiftDict[1] = zero(GCC.H2)
  #first, assume basePoint = 1
  while length(keys(shiftDict)) < nv # this is why we checked connectedness of G first.
    for v1 in keys(shiftDict)
      for v2 in all_neighbors(G.g, v1)
        haskey(shiftDict, v2) && continue
        e = Edge(v1,v2)
        we = G.w[e]
        wz = sum([we[i] * weight[i] for i in 1:r])
        shiftDict[v2] = shiftDict[v1] -wz * edgeCurveClass(G, e)
      end
    end
  end
  # Now shift so that the correct base point has zero vertical class.
  horToGH2 = ModuleHomomorphism(horizontalEdgeLattice, GCC.H2, [shiftDict[i] - shiftDict[basePoint] for i in 1:nv])
  qMatrix = [Gq Int(0); Gq Int(0); horToGH2 horToExtra ]
  Sq = ModuleHomomorphism(SedgeLattice, SH2, qMatrix)
  SedgeToGenIndex = Dict{Edge, Int64}()
  for e in edges(G.g)
    indE = GCC.edgeToGenIndex[e]
    SedgeToGenIndex[e] = indE
    SedgeToGenIndex[reverse(e)] = indE
    SedgeToGenIndex[Edge(src(e)+nv, dst(e)+nv)] = indE + ne
    SedgeToGenIndex[Edge(dst(e)+nv, src(e)+nv)] = indE + ne
  end
  for v in 1:nv
    SedgeToGenIndex[Edge(v, v+nv)] = 2*ne + v
    SedgeToGenIndex[Edge(v+nv, v)] = 2*ne + v
  end
  dualConeRaySum, C, H2ToCN = _finish_GKM_H2(SedgeLattice, SH2, Sq, SG, SedgeToGenIndex)
  SG.curveClasses = GKM_H2(SG, SedgeLattice, SH2, SedgeToGenIndex, Sq, dualConeRaySum, C, H2ToCN, nothing, projs[1])

  return SG
end

function _SeidelSectionCount(SG::AbstractGKM_graph)
  @req divides(n_vertices(SG.g), 2)[1] "SG is not a Seidel space!"
  H2 = GKM_second_homology(SG)
  rkH2 = length(gens(H2.H2)) - 1

  if !isnothing(H2.sectionCount)
    return H2.sectionCount
  end

  # generate sectional multiplicity as ZZ-module homorphism from H2 to ZZ
  ZZasModule = free_module(ZZ, 1)
  H2.sectionCount = ModuleHomomorphism(H2.H2, ZZasModule, vcat([zero(ZZasModule) for i in 1:rkH2], [gens(ZZasModule)[1]]))
  return H2.sectionCount
end

"""
Warning: Don't use this on anything that is not the output of Seidel_space(...).
"""
function _effectiveSectionClassesWithChernNumber(
  SG::AbstractGKM_graph,
  chernNumber::ZZRingElem;
)

  H2 = GKM_second_homology(SG)
  secCt = _SeidelSectionCount(SG)
  cN = H2.chernNumber

  ZZ2, _, _ = direct_sum([codomain(cN), codomain(secCt)])
  q = ModuleHomomorphism(H2.H2, ZZ2, hcat(matrix(cN), matrix(secCt)))
  e0 = nothing
  try
    e0 = preimage(q, chernNumber * gens(ZZ2)[1] + gens(ZZ2)[2])
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