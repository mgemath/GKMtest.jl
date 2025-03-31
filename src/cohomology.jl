"""
    _equivariant_cohomology_ring(G::AbstractGKM_graph)

Return the equivariant cohomology ring of the GKM graph. This function is called internally once during
construction of the GKm graph.
"""
function _equivariant_cohomology_ring(G::AbstractGKM_graph)::GKM_cohomology_ring
  coeffRing, _ = polynomial_ring(QQ, ["t$i" for i in 1:rank_torus(G)])
  coeffRingLocalized = fraction_field(coeffRing)
  cohomRing = free_module(coeffRing, n_vertices(G.g))
  cohomRingLocalized = free_module(coeffRingLocalized, n_vertices(G.g))
  edgeWeightClasses = Dict{Edge, QQMPolyRingElem}()
  pointEulerClasses = vcat(Union{Nothing, QQMPolyRingElem}[], repeat([nothing], n_vertices(G.g)))

  return GKM_cohomology_ring(G, coeffRing, coeffRingLocalized, cohomRing, cohomRingLocalized, edgeWeightClasses, pointEulerClasses)
end

@doc raw"""
    is_gkm_class(c::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph) -> Bool

Return true if the given class is an actual GKM cohomology class.
This holds if and only if the difference between localizations at fixed points connected through
an edge e is divisible by the weight of e.
"""
function is_gkm_class(c::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph)::Bool
  return is_gkm_class(c, G.equivariantCohomology)
end

@doc raw"""
    is_gkm_class(c::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring) -> Bool

Return true if the given class is an actual GKM cohomology class.
This holds if and only if the difference between localizations at fixed points connected through
an edge e is divisible by the weight of e.
"""
function is_gkm_class(c::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::Bool

  g = R.gkm.g
  w = R.gkm.w

  for e in edges(g)
    polySrc = c[src(e)]
    polyDst = c[dst(e)]
    dif = polySrc - polyDst
    if !divides(dif, weight_class(e, R))[1]
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

@doc raw"""
    point_class(vertexLabel::String, G::AbstractGKM_graph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the fixed point with given label
"""
function point_class(vertexLabel::String, G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}
  return point_class(vertexLabel, G.equivariantCohomology)
end

function point_class(G::AbstractGKM_graph, vertexLabel::String)::FreeModElem{QQMPolyRingElem}
  return point_class(vertexLabel, G.equivariantCohomology)
end

function point_class(vertexLabel::String, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  @req vertexLabel in R.gkm.labels "Vertex not found"
  return point_class(indexin([vertexLabel], R.gkm.labels)[1], R)
end

@doc raw"""
    point_class(vertex::Int, G::AbstractGKM_graph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the given fixed point
"""
function point_class(vertex::Int, G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}
  return point_class(vertex, G.equivariantCohomology)
end 

function point_class(G::AbstractGKM_graph, vertex::Int)::FreeModElem{QQMPolyRingElem}
  return point_class(vertex, G.equivariantCohomology)
end 

function point_class(vertex::Int, R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  return euler_class(vertex, R) * gens(R.cohomRing)[vertex]
end 

@doc raw"""
    poincare_dual(gkmSub::AbstractGKM_subgraph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual cohomology class of the GKM subgraph.
"""
function poincare_dual(gkmSub::AbstractGKM_subgraph)::FreeModElem{QQMPolyRingElem}

  R = gkmSub.super.equivariantCohomology

  res = zero(R.cohomRing)

  for v in gkmSub.vDict
    vContrib = R.coeffRing(1)
    for w in all_neighbors(gkmSub.super.g, v)
      e = Edge(v,w)
      if !has_edge(gkmSub, e)
        vContrib *= weight_class(e, R)
      end
    end
    res += vContrib * gens(R.cohomRing)[v]
  end

  return res
end

@doc raw"""
    weight_class(e::Edge, G::AbstractGKM_graph) -> QQMPolyRingElem

Return the weight of the edge e as an element of the coefficient ring of the equivariant cohomology theory.
"""
function weight_class(e::Edge, G::AbstractGKM_graph)::QQMPolyRingElem
  return weight_class(e, G.equivariantCohomology)
end

@doc raw"""
    weight_class(e::Edge, R::GKM_cohomology_ring) -> QQMPolyRingElem

Return the weight of the edge e as an element of the coefficient ring of the equivariant cohomology theory.
"""
function weight_class(e::Edge, R::GKM_cohomology_ring)::QQMPolyRingElem
  try
    return R.edgeWeightClasses[e]
  catch err
    if isa(err, KeyError)
      res = _weight_class(e, R)
      R.edgeWeightClasses[e] = res
      return res
    else
      rethrow(err)
    end
  end
end

function _weight_class(e::Edge, R::GKM_cohomology_ring)::QQMPolyRingElem
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

@doc raw"""
    euler_class(vertex::Int, G::AbstractGKM_graph) -> QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point v.
"""
function euler_class(vertex::Int, G::AbstractGKM_graph)::QQMPolyRingElem
  return euler_class(vertex, G.equivariantCohomology)
end

@doc raw"""
    euler_class(vertex::Int, R::GKM_cohomology_ring) -> QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point v.
"""
function euler_class(vertex::Int, R::GKM_cohomology_ring)::QQMPolyRingElem
  res = R.pointEulerClasses[vertex]
  if isnothing(res)
    res = _euler_class(vertex, R)
    R.pointEulerClasses[vertex] = res
  end
  return res
end

function _euler_class(vertex::Int, R::GKM_cohomology_ring)::QQMPolyRingElem
  res = R.coeffRing(1)
  for i in all_neighbors(R.gkm.g, vertex)
    res = mul!(res, weight_class(Edge(vertex, i), R))
  end
  return res
end

@doc raw"""
    euler_class(vertexLabel::String, G::AbstractGKM_graph) -> QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point with the given label.
"""
function euler_class(vertexLabel::String, G::AbstractGKM_graph)::QQMPolyRingElem
  return euler_class(vertexLabel, G.equivariantCohomology)
end

@doc raw"""
    euler_class(vertexLabel::String, R::GKM_cohomology_ring) -> QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point with the given label.
"""
function euler_class(vertexLabel::String, R::GKM_cohomology_ring)::QQMPolyRingElem
  @req vertexLabel in R.gkm.labels "Vertex not found"
  return euler_class(indexin([vertexLabel], R.gkm.labels)[1], R)
end

@doc raw"""
    integrate_gkm_class(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph; check::Bool=true) -> QQMPolyRingElem

Integrate the GKM class, yielding an element of the coefficient ring. This only works if `is_gkm_class(class,R) == true`.
"""
function integrate_gkm_class(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph; check::Bool=true)::QQMPolyRingElem
  return integrate_gkm_class(class, G; check)
end

@doc raw"""
    integrate_gkm_class(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring; check::Bool=true) -> QQMPolyRingElem

Integrate the GKM class, yielding an element of the coefficient ring. This only works if `is_gkm_class(class,R) == true`.
"""
function integrate_gkm_class(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring; check::Bool=true)::QQMPolyRingElem
  if check
    @req is_gkm_class(class, R) "Given class is not GKM."
  end
  res = integrate(class, R)
  (flag, quotient) = divides(numerator(res), denominator(res))
  if check
    @req flag "Integral of GKMclass is fraction"
  end
  return quotient
end


@doc raw"""
    first_chern_class(G::AbstractGKM_graph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant first chern class of the GKM space.
"""
function first_chern_class(G::AbstractGKM_graph)::FreeModElem{QQMPolyRingElem}
  return first_chern_class(G.equivariantCohomology)
end

@doc raw"""
    first_chern_class(R::GKM_cohomology_ring) -> FreeModElem{QQMPolyRingElem}

Return the equivariant first chern class of the GKM space of which R is the cohomology ring.
"""
function first_chern_class(R::GKM_cohomology_ring)::FreeModElem{QQMPolyRingElem}
  res = zero(R)
  for v in 1:n_vertices(R.gkm.g)
    localFactor = zero(R.coeffRing)
    for w in all_neighbors(R.gkm.g, v)
      localFactor += weight_class(Edge(v, w), R)
    end
    res += localFactor * gens(R.cohomRing)[v]
  end
  return res
end

@doc raw"""
    integrate(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph, e::Edge) -> AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class over the curve represented by the GKM graph edge.
"""
function integrate(
  class::FreeModElem{QQMPolyRingElem},
  G::AbstractGKM_graph,
  e::Edge
)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  return integrate(class, G.equivariantCohomology, e)
end

@doc raw"""
    integrate(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring, e::Edge) -> AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class over the curve represented by the GKM graph edge.
"""
function integrate(
  class::FreeModElem{QQMPolyRingElem},
  R::GKM_cohomology_ring,
  e::Edge
)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  return (class[src(e)] - class[dst(e)]) // weight_class(e, R)
end

@doc raw"""
    integrate(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph) -> AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class, yielding an element of the fraction field of the equivariant coefficient ring.
If `is_gkm_class(class, R) == true`, this fraction will represent a ring element.
Use `integrate_gkm_class()` instead if you know that the class is a GKM class.
"""
function integrate(class::FreeModElem{QQMPolyRingElem}, G::AbstractGKM_graph)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  return integrate(class, G.equivariantCohomology)
end

@doc raw"""
    integrate(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring) -> AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class, yielding an element of the fraction field of the equivariant coefficient ring.
If `is_gkm_class(class, R) == true`, this fraction will represent a ring element.
Use `integrate_gkm_class()` instead if you know that the class is a GKM class.
"""
function integrate(class::FreeModElem{QQMPolyRingElem}, R::GKM_cohomology_ring)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  res = R.coeffRing(0)
  for v in 1:n_vertices(R.gkm.g)
    euler = euler_class(v, R)
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