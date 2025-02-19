
import Base: zero, one, *, ^

"""
    _equivariant_cohomology_ring(G::AbstractGKM_graph)

Return the equivariant cohomology ring of the GKM graph. This function is called internally once during
construction of the GKm graph.
"""
function _equivariant_cohomology_ring(G::AbstractGKM_graph)::GKM_cohomology_ring
  coeffRing, _ = polynomial_ring(QQ, ["t$i" for i in 1:rank_torus(G)])
  cohomRing = free_module(coeffRing, n_vertices(G.g))
  edgeWeightClasses = Dict{Edge, QQMPolyRingElem}()
  pointEulerClasses = QQMPolyRingElem[]
  sizehint!(pointEulerClasses, n_vertices(G.g))

  res = GKM_cohomology_ring(G, coeffRing, cohomRing, edgeWeightClasses, pointEulerClasses)

  for e in edges(G.g)
    res.edgeWeightClasses[e] = _weightClass(e, res)
    res.edgeWeightClasses[reverse(e)] = _weightClass(reverse(e), res)
  end
  for v in 1:n_vertices(G.g)
    push!(res.pointEulerClasses, _eulerClass(v, res))
  end
  return res
end

"""
    isGKMclass(c::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph)::Bool

Return true if the given class is an actual GKM cohomology class.
This holds if and only if the difference between localizations at fixed points connected through
an edge e is divisible by the weight of e.
"""
function isGKMclass(c::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph)::Bool
  return isGKMclass(c, G.equivariantCohomology)
end

"""
    isGKMclass(c::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::Bool

Return true if the given class is an actual GKM cohomology class.
This holds if and only if the difference between localizations at fixed points connected through
an edge e is divisible by the weight of e.
"""
function isGKMclass(c::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::Bool

  g = R.gkm.g
  w = R.gkm.w

  for e in edges(g)
    polySrc = c[src(e)]
    polyDst = c[dst(e)]
    dif = polySrc - polyDst
    if !divides(dif, weightClass(e, R))[1]
      return false
    end
  end
  return true
end

function scalar(s::Union{Int, Rational, QQFieldElem}, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  z = R.coeffRing(s)
  return R.cohomRing([z for _ in 1:n_vertices(R.gkm.g)])
end

function scalar(s::QQMPolyRingElem, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  return R.cohomRing([s for _ in 1:n_vertices(R.gkm.g)])
end 

function zero(R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  return scalar(0,R)
end

function one(R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  return scalar(1,R)
end

"""
    pointClass(vertexLabel::String, G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the fixed point with given label
"""
function pointClass(vertexLabel::String, G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}
  return pointClass(vertexLabel, G.equivariantCohomology)
end

"""
    pointClass(vertexLabel::String, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the fixed point with given label
"""
function pointClass(vertexLabel::String, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  @req vertexLabel in R.gkm.labels "Vertex not found"
  return pointClass(indexin([vertexLabel], R.gkm.labels)[1], R)
end

"""
    pointClass(vertex::Int, G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the given fixed point
"""
function pointClass(vertex::Int, G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}
  return pointClass(vertex, G.equivariantCohomology)
end 

"""
    pointClass(vertex::Int, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the given fixed point
"""
function pointClass(vertex::Int, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  return eulerClass(vertex, R) * gens(R.cohomRing)[vertex]
end 

"""
    PDClass(gkmSub::AbstractGKM_subgraph)::FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual cohomology class of the GKM subgraph.
"""
function PDClass(gkmSub::AbstractGKM_subgraph)::FreeModElem{QQMPolyRingElem}

  R = gkmSub.super.equivariantCohomology

  res = zero(R.cohomRing)

  for v in gkmSub.vDict
    vContrib = R.coeffRing(1)
    for w in all_neighbors(gkmSub.super.g, v)
      e = Edge(v,w)
      if !has_edge(gkmSub, e)
        vContrib *= weightClass(e, R)
      end
    end
    res += vContrib * gens(R.cohomRing)[v]
  end

  return res
end

"""
    weightClass(e::Edge, G::AbstractGKM_graph)::QQMPolyRingElem

Return the weight of the edge e as an element of the coefficient ring of the equivariant cohomology theory.
"""
function weightClass(e::Edge, G::AbstractGKM_graph)::QQMPolyRingElem
  return weightClass(e, G.equivariantCohomology)
end

"""
    weightClass(e::Edge, R::GKM_cohomology_ring)::QQMPolyRingElem

Return the weight of the edge e as an element of the coefficient ring of the equivariant cohomology theory.
"""
function weightClass(e::Edge, R::GKM_cohomology_ring)::QQMPolyRingElem
  return R.edgeWeightClasses[e]
end

function _weightClass(e::Edge, R::GKM_cohomology_ring)::QQMPolyRingElem
  w = R.gkm.w[e]
  rk = rank_torus(R.gkm)
  coeffs = R.coeffRing
  t = gens(coeffs)
  
  res = coeffs(0)
  for i in 1:rk
    res += w[i] * t[i]
  end
  return res
end

"""
    eulerClass(vertex::Int, G::AbstractGKM_graph)::QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point v.
"""
function eulerClass(vertex::Int, G::AbstractGKM_graph)::QQMPolyRingElem
  return eulerClass(vertex, G.equivariantCohomology)
end

"""
    eulerClass(vertex::Int, R::GKM_cohomology_ring)::QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point v.
"""
function eulerClass(vertex::Int, R::GKM_cohomology_ring)::QQMPolyRingElem
  return R.pointEulerClasses[vertex]
end

function _eulerClass(vertex::Int, R::GKM_cohomology_ring)::QQMPolyRingElem
  res = R.coeffRing(1)
  for i in all_neighbors(R.gkm.g, vertex)
    res = mul!(res, weightClass(Edge(vertex, i), R))
  end
  return res
end

"""
    eulerClass(vertexLabel::String, G::AbstractGKM_graph)::QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point with the given label.
"""
function eulerClass(vertexLabel::String, G::AbstractGKM_graph)::QQMPolyRingElem
  return eulerClass(vertexLabel, G.equivariantCohomology)
end

"""
    eulerClass(vertexLabel::String, R::GKM_cohomology_ring)::QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point with the given label.
"""
function eulerClass(vertexLabel::String, R::GKM_cohomology_ring)::QQMPolyRingElem
  @req vertexLabel in R.gkm.labels "Vertex not found"
  return eulerClass(indexin([vertexLabel], R.gkm.labels)[1], R)
end

"""
    integrateGKMclass(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph; check::Bool=true)::QQMPolyRingElem

Integrate the GKM class, yielding an element of the coefficient ring. This only works if isGKMclass(class,R) == true.
"""
function integrateGKMClass(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph; check::Bool=true)::QQMPolyRingElem
  return integrateGKMClass(class, G; check)
end

"""
    integrateGKMclass(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring; check::Bool=true)::QQMPolyRingElem

Integrate the GKM class, yielding an element of the coefficient ring. This only works if isGKMclass(class,R) == true.
"""
function integrateGKMClass(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring; check::Bool=true)::QQMPolyRingElem
  if check
    @req isGKMclass(class, R) "Given class is not GKM."
  end
  res = integrateClass(class, R)
  (flag, quotient) = divides(numerator(res), denominator(res))
  if check
    @req flag "Integral of GKMclass is fraction"
  end
  return quotient
end


"""
Return the equivariant first chern class of the GKM space.
"""
function firstChernClass(G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}
  return firstChernClass(G.equivariantCohomology)
end

"""
Return the equivariant first chern class of the GKM space of which R is the cohomology ring.
"""
function firstChernClass(R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  res = zero(R)
  for v in 1:n_vertices(R.gkm.g)
    localFactor = zero(R.coeffRing)
    for w in all_neighbors(R.gkm.g, v)
      localFactor += weightClass(Edge(v, w), R)
    end
    res += localFactor * gens(R.cohomRing)[v]
  end
  return res
end

"""
Integrate the cohomology class over the curve represented by the GKM graph edge.
"""
function integrateOverEdge(
  class::FreeModElem{QQMPolyRingElem},
  G::AbstractGKM_graph,
  e::Edge
)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  return integrateOverEdge(class, G.equivariantCohomology, e)
end

"""
Integrate the cohomology class over the curve represented by the GKM graph edge.
"""
function integrateOverEdge(
  class::FreeModElem{QQMPolyRingElem},
  R::GKM_cohomology_ring,
  e::Edge
)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  return (class[src(e)] - class[dst(e)]) // weightClass(e, R)
end

"""
    integrateClass(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class, yielding an element of the fraction field of the equivariant coefficient ring.
If isGKMclass(class, R) == true, this fraction will represent a ring element.
Use integrateGKMClass() instead if you know that the class is a GKM class.
"""
function integrateClass(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  return integrateClass(class, G.equivariantCohomology)
end

"""
    integrateClass(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class, yielding an element of the fraction field of the equivariant coefficient ring.
If isGKMclass(class, R) == true, this fraction will represent a ring element.
Use integrateGKMClass() instead if you know that the class is a GKM class.
"""
function integrateClass(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  res = R.coeffRing(0)
  for v in 1:n_vertices(R.gkm.g)
    euler = eulerClass(v, R)
    contrib = class[v] // euler
    res = add!(res, contrib)
  end
  return res
end

function multiply(a::FreeModElem{QQMPolyRingElem}, b::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  return R.cohomRing([a[i] * b[i] for i in 1:n_vertices(R.gkm.g)])
end

function *(a::FreeModElem{T}, b::FreeModElem{T})::FreeModElem{T} where T <: RingElem
  return parent(a)([a[i] * b[i] for i in 1:rank(parent(a))])
end

function ^(a::FreeModElem{T}, n::Int64)::FreeModElem{T} where T <: RingElem
  return parent(a)([a[i]^n for i in 1:rank(parent(a))])
end

function Base.show(io::IO, R::GKM_cohomology_ring)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM cohomology ring")
  else
    # nested printing allowed, preferably terse
    print(io, "Cohomology ring of GKM graph with $(n_vertices(R.gkm.g)) nodes, valency $(valency(R.gkm)), and rank $(rank_torus(R.gkm)) torus")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", R::GKM_cohomology_ring)

  print(io, "Cohomology ring of GKM graph with $(n_vertices(R.gkm.g)) nodes, valency $(valency(R.gkm)), and rank $(rank_torus(R.gkm)) torus\n")
  print(io, "Coefficient Ring: ")
  show(io, MIME"text/plain"(), R.coeffRing)
  print(io, "\nCohomology Ring: ")
  show(io, MIME"text/plain"(), R.cohomRing)
  
end