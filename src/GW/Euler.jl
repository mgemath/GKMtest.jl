# Warning: This is not actually the inverse of the Euler class, as the h classes will be multiplied later.
function Euler_inv(dt::GW_decorated_tree)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  C = dt.gkm.equivariantCohomology.coeffRing
  res = C(1)//C(1)

  for v in 1:n_vertices(dt.tree)

    valv = degree(dt.tree, v)
    res = res * (eulerClass(imageOf(v, dt), dt.gkm))^(valv - 1)

    tmpSum = C(0)//C(1)

    for v2 in all_neighbors(dt.tree, v)
      e = Edge(v,v2)
      wev = weightClass(imageOf(e, dt), dt.gkm) // edgeMult(e, dt)
      res = res // wev
      tmpSum = tmpSum + 1//wev
    end

    res = res * tmpSum^( valv - 3 +  count(i -> i==v, dt.marks) )
  end

  return res
end

"""
Calculate h(epsilon, d) as in [Liu--Sheshmani, Lemma 4.5, p. 16].
"""
function _h(e::Edge, d::Int, con::GKM_connection, R::GKM_cohomology_ring; check::Bool=true)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  gkm = con.gkm
  C = R.coeffRing

  if check
    @req con.gkm == R.gkm "GKM connection and cohomology ring don't belong to the same GKM graph"
    @req has_edge(gkm.g, e) "edge not found in GKM graph"
    @req d>0 "d is non-positive"
  end

  we = weightClass(e, R) # weight of the edge e

  res = ( C((-1)^d) * ZZ(d)^(2d) ) // ( factorial(ZZ(d))^2 )
  res = res // ( (we)^(2d) )

  for v in all_neighbors(gkm.g, src(e))

    v == dst(e) && continue 

    ei = Edge(src(e), v)
    wei = weightClass(ei, R)
    ai = con.a[(e, ei)]
    res = res * _b(1//d * we, wei, d*ai, C)
  end

  return res
end

"""
Calculate b(u,w,a) as in [Liu--Sheshmani, Lemma 4.5, p.16]. C is the coefficient ring
"""
function _b(u::QQMPolyRingElem, w::QQMPolyRingElem, a::ZZRingElem, C::QQMPolyRing)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  res = C(1) // C(1) # make sure this has FracFieldElem type.
  if a >= 0
    for j in 0:a
      res = res // (w - j*u)
    end
  else 
    for j in 1:(-a-1)
      res = res // (w + j*u)
    end
  end
  return res
end

function GWTreeContribution(
  dt::GW_decorated_tree,
  P_input;
  check::Bool=true)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  res = Euler_inv(dt)
  R = dt.gkm.equivariantCohomology
  con = get_GKM_connection(dt.gkm)

  # multiply by h classes
  for e in edges(dt.tree)
    res *= _h(imageOf(e, dt), edgeMult(e, dt), con, R; check)
  end

  #multiply by input class
  res *= Base.invokelatest(P_input.func, dt)

  return res
end