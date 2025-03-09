export flag_variety, grassmannian, gkm_graph_of_toric
# import Oscar: projective_space

@doc raw"""
    flag_variety(::Type{GKM_graph}, s::Vector{Int64}) -> AbstractGKM_graph{ZZRingElem}

Construct the GKM graph of the variety of flags of ``\mathbb{C}^n``. The dimensions of quotients are expressed by the array `s`.

# Examples
```jldoctest
julia> flag_variety(GKM_graph, [2,1])
GKM graph with 3 nodes, valency 2 and axial function:
13 -> 12 => (0, -1, 1)
23 -> 12 => (-1, 0, 1)
23 -> 13 => (-1, 1, 0)

```
"""
function flag_variety(::Type{GKM_graph}, s::Vector{Int64})

    @req !isempty(s) "the vector of dimensions is empty"
    @req all(i -> s[i] > 0, eachindex(s)) "all dimensions must be positive"
  
    K = [sum(s[1:i]) for i in 0:length(s)]
    d = Dict{Int64, Vector{Int64}}()
    index = 1
    
    for c in Combinatorics.permutations(1:K[end])
      if all(i-> issorted(c[(K[i]+1):K[i+1]]), 1:(length(K)-1))
        d[index] = c
        index += 1
      end
    end
        
    nv = length(keys(d))
    g = Graph{Undirected}(nv)
  
    M = free_module(ZZ, K[end])
    W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
      
    for v in 1:nv
      for w in (v+1):nv
        dif = _T_link(d[v], d[w], K)
          
        isempty(dif) && continue
          
        add_edge!(g, v, w)
        W[Edge(w, v)] = gens(M)[dif[2]] - gens(M)[dif[1]]
      end
    end
  
    labels = [prod(i -> d[n][i]>9 ? "|$(d[n][i])|" : "$(d[n][i])", 1:K[end-1]) for n in 1:nv] 

    return gkm_graph(g, labels, M, W)
  end
  
  
# this function variefies if two flags can be connected by a T-invariant curve. If this is not possible, an empty array is returned. Otherwise, it returns the two indices that need to be swapped
function _T_link(a::Vector{Int64}, b::Vector{Int64}, K::Vector{Int64})
    
  ans = Int64[]
  
  for i in 1:length(K)
    a1 = a[1:K[i]]
    b1 = b[1:K[i]]
    a1_b1 = setdiff(a1, b1)

    if isempty(a1_b1)
      continue
    end

    if length(a1_b1) > 1
      return Int64[]
    end

    b1_a1 = setdiff(b1, a1)

    if isempty(ans)
      ans = [a1_b1[1], b1_a1[1]]
      continue
    end

    if ans != [a1_b1[1], b1_a1[1]]
      return Int64[]
    end
end

  return ans
end

@doc raw"""
    grassmannian(::Type{gkm_graph}, k::Int, n::Int) -> AbstractGKM_graph{ZZRingElem}

Construct the Grassmann variety of `k`-planes in the complex vector space of dimension `n`.

# Examples
```jldoctest
julia> grassmannian(GKM_graph, 1,4)
GKM graph with 5 nodes, valency 4 and axial function:
2 -> 1 => (-1, 1, 0, 0, 0)
3 -> 1 => (-1, 0, 1, 0, 0)
3 -> 2 => (0, -1, 1, 0, 0)
4 -> 1 => (-1, 0, 0, 1, 0)
4 -> 2 => (0, -1, 0, 1, 0)
4 -> 3 => (0, 0, -1, 1, 0)
5 -> 1 => (-1, 0, 0, 0, 1)
5 -> 2 => (0, -1, 0, 0, 1)
5 -> 3 => (0, 0, -1, 0, 1)
5 -> 4 => (0, 0, 0, -1, 1)

```
"""
function grassmannian(::Type{GKM_graph}, k::Int, n::Int)
  @req (k >= 0 && n >= k) "Dimension must be non-negative"
  
  return flag_variety(GKM_graph, [k, n])
end

# @doc raw"""
#     projective_space(::Type{gkm_graph}, d::Int)

# Construct the projective space of dimension `d`.

# # Examples
# ```jldoctest
# julia> projective_space(NormalToricVariety, 2)
# Normal toric variety
# ```
# """
# function projective_space(::Type{GKM_graph}, d::Int)
#   @req d >= 0 "Dimension must be non-negative"
  
#   return grassmannian(GKM_graph, 1, d+1)
# end


@doc raw"""
    gkm_graph_of_toric(v::NormalToricVariety) -> AbstractGKM_graph{ZZRingElem}

Construct the GKM graph of the (smooth, projective) toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> gkm_graph_of_toric(P2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1)
3 -> 1 => (0, 1, -1)
3 -> 2 => (-1, 1, 0)

julia> F = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> gkm_graph_of_toric(F)
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1, 0)
3 -> 2 => (3, 1, 0, -1)
4 -> 1 => (0, 1, 3, -1)
4 -> 3 => (-1, 0, 1, 0)
```
"""
function gkm_graph_of_toric(v::NormalToricVariety)

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
          ud = lcm(denominator.(pol_ray)) * pol_ray
          break
      end
      break
  end

  ans = zero(M)

  for (k, vi) in enumerate(rays(v))
      ans += Int64(dot(vi, ud))*scalars[k]
  end

  return ans
end