#TODO: If calculating GKM_second_homology of the Seidel space takes too long, implement it here.

function Seidel_space(
  G::AbstractGKM_graph,
  w::AbstractAlgebra.Generic.FreeModuleElem{R}
) where R<:GKM_weight_type

  @req parent(w) === G.M "The weight does not belong to the right character lattice"

  nv = n_vertices(G.g)
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
    GKMadd_edge!(SG, i, i+nv, z)
  end
  for e in edges(G.g)
    # edge over 0:
    GKMadd_edge!(SG, src(e), dst(e), iM(G.w[e]))
    # edge over inf:
    wz = sum([G.w[e][i] * w[i] for i in 1:r])
    GKMadd_edge!(SG, nv + src(e), nv + dst(e), iM(G.w[e]) -  wz * z)
  end

  # infer GKM connection to Seidel space if possible
  con = get_GKM_connection(G)
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
    set_GKM_connection!(SG, Scon)
  end

  return SG
end

function _SeidelSectionCount(SG::AbstractGKM_graph)
  @req divides(n_vertices(SG.g), 2)[1] "SG is not a Seidel space!"
  nv = div(n_vertices(SG.g), 2)
  H2 = GKM_second_homology(SG)

  if !isnothing(H2.sectionCount)
    return H2.sectionCount
  end

  # generate sectional multiplicity as ZZ-module homorphism from H2 to ZZ
  ZZasModule = free_module(ZZ, 1)
  edgeSectionDict = Dict{Int64, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  edgeSectionVect = Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for e in edges(SG.g)
    if abs(src(e) - dst(e)) == nv
      edgeSectionDict[H2.edgeToGenIndex[e]] = gens(ZZasModule)[1]
    else
      edgeSectionDict[H2.edgeToGenIndex[e]] = 0 * gens(ZZasModule)[1]
    end
  end
  for i in 1:n_edges(SG.g)
    push!(edgeSectionVect, edgeSectionDict[i])
  end
  edgeLatticeToSec = ModuleHomomorphism(H2.edgeLattice, ZZasModule, edgeSectionVect)

  H2ToSec = Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for g in gens(H2.H2)
    p = preimage(H2.quotientMap, g)
    push!(H2ToSec, edgeLatticeToSec(p))
  end
  H2.sectionCount = ModuleHomomorphism(H2.H2, ZZasModule, H2ToSec)
  return H2.sectionCount
end

"""
Warning: Don't use this on anything that is not the output of Seidel_space(...).
"""
function _effectiveSectionClassesWithChernNumber(
  SG::AbstractGKM_graph,
  chernNumber::ZZRingElem;
)

  # TODO: this whole thing!
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