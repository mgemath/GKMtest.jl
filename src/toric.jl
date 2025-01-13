#this function serves to convert a toric variety into a GKM variety

function gkm_graph(v::NormalToricVariety)

  @req is_smooth(v) && is_projective(v) "toric variety must be smooth and projective"

  len = length(maximal_cones(v))
  g = Graph{Undirected}(len)
  M = free_module(ZZ, n_rays(v))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()

  for sigma1 in 1:(len-1)
    for sigma2 in (sigma1+1):len
      count(x -> x in rays(maximal_cones(v)[sigma1]), rays(maximal_cones(v)[sigma2])) != (dim(v) - 1) && continue
      add_edge!(g, sigma1, sigma2)
      W[Edge(sigma2, sigma1)] = _omega(v, sigma1, sigma2, M)
      # W[Edge(sigma2, sigma1)] = -W[Edge(sigma1, sigma2)]
    end
  end

  return gkm_graph(g, ["$i" for i in 1:len], M, W)
end

function _omega(v::NormalToricVariety, n_SIGMA1::Int64, n_SIGMA2::Int64, M::AbstractAlgebra.Generic.FreeModule{ZZRingElem})::AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}

    SIGMA1 = maximal_cones(v)[n_SIGMA1]
    SIGMA2 = maximal_cones(v)[n_SIGMA2]
    scalars = gens(M)
    ud = QQFieldElem[]

    for ray in rays(SIGMA1)
        ray in rays(SIGMA2) && continue
        for pol_ray in rays(polarize(SIGMA1))
            dot(ray, pol_ray) == 0 && continue
            # ud = pol_ray
            ud = lcm(denominator.(pol_ray)) * pol_ray
            break
        end
        break
    end

    # ans = F(0)
    ans = zero(M)

    for (k, vi) in enumerate(rays(v))# keys(rays_enum)
        ans += Int64(dot(vi, ud))*scalars[k]
        # ans += numerator(dot(vi, ud)) * scalars[k]
    end

    return ans
end