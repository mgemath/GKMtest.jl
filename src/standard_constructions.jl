export flag_variety, grassmannian, gkm_graph_of_toric, projective_space, schubert_class, schubert_classes

@doc raw"""
    flag_variety(::Type{GKM_graph}, s::Vector{Int64}) -> AbstractGKM_graph{ZZRingElem}

Construct the GKM graph of the variety of flags of ``\mathbb{C}^n``. The dimensions of quotients are expressed by the array `s`. The labels represent the vectors generating the flags. For example, if ``s=[1,2,1]``, the string ``213`` corresponds to the flag:

``0\subset \langle e_2 \rangle \subset \langle e_2, e_1, e_3 \rangle \subset \langle e_2, e_1, e_3, e_4 \rangle=\mathbb{C}^4.``

!!! note
    This function is faster than `generalized_gkm_flag(root_system(:A, n-1), S)`, but the results are isomorphic.

# Examples
```jldoctest
julia> flag_variety(GKM_graph, [1,3])
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

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

  K::Vector{Int64} = [sum(s[1:i]) for i in 0:length(s)]
  d::Dict{Int64, NTuple{K[end], Int64}} = Dict{Int64, NTuple{K[end], Int64}}()
  index::Int64 = 1
  
  for c in Combinatorics.permutations(1:K[end])
    if all(i-> issorted(c[(K[i]+1):K[i+1]]), 1:length(s))
      d[index] = (c...,)
      index += 1
    end
  end
      
  nv::Int64 = length(keys(d))
  g = Graph{Undirected}(nv)

  M = free_module(ZZ, K[end])
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
    
  get_root::Dict{Edge, Vector{Int64}} = Dict{Edge, Vector{Int64}}()

  for v in 1:nv
    for w in (v+1):nv
      dif = _T_link(d[v], d[w], K)
        
      isempty(dif) && continue

      add_edge!(g, v, w)
      W[Edge(w, v)] = gens(M)[dif[2]] - gens(M)[dif[1]]

      get_root[Edge(w, v)] = dif
    end
  end

  ## construct connection
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} = Dict{Tuple{Edge, Edge}, ZZRingElem}()

  for _v in vertices(g)
    for _w in all_neighbors(g, _v)

      alpha = get_root[_w < _v ? Edge(_v, _w) : Edge(_w, _v)]


      for _u in all_neighbors(g, _v)
        # e = (_v, _w)
        # e'= (_v, _u)

        beta = get_root[_u < _v ? Edge(_v, _u) : Edge(_u, _v)]

        if alpha == beta
          a[(Edge(_v, _w), Edge(_v, _u))] = ZZ(2)
        elseif alpha[1] == beta[1] || alpha[2] == beta[2]
          a[(Edge(_v, _w), Edge(_v, _u))] = ZZ((_v < _w ? -1 : 1)*(_v < _u ? -1 : 1))
        elseif alpha[1] == beta[2] || alpha[2] == beta[1]
          a[(Edge(_v, _w), Edge(_v, _u))] = -ZZ((_v < _w ? -1 : 1)*(_v < _u ? -1 : 1))
        else
          a[(Edge(_v, _w), Edge(_v, _u))] = ZZ(0)
        end
        
      end
    end
  end

  labels = [prod(i -> d[n][i]>9 ? "|$(d[n][i])|" : "$(d[n][i])", 1:K[end-1]) for n in 1:nv] 

  G = gkm_graph(g, labels, M, W)


  con = build_GKM_connection(G, a)
  set_connection!(G, con)

  return G

end

  
# this function variefies if two flags can be connected by a T-invariant curve. If this is not possible, an empty array is returned. Otherwise, it returns the two indices that need to be swapped
function _T_link(a::Tuple{Vararg{Int64}}, b::Tuple{Vararg{Int64}}, K::Vector{Int64})
    
  ans::Vector{Int64} = [0, 0]
  
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

    if ans[1] == 0
      ans = [a1_b1[1], b1_a1[1]]
      continue
    end

    if ans != [a1_b1[1], b1_a1[1]]
      return Int64[]
    end
end

  return sort(ans)
end

@doc raw"""
    grassmannian(::Type{gkm_graph}, k::Int, n::Int) -> AbstractGKM_graph{ZZRingElem}

Construct the Grassmann variety of `k`-planes in the complex vector space of dimension `n`.

# Examples
```jldoctest
julia> grassmannian(GKM_graph, 2, 4)
GKM graph with 6 nodes, valency 4 and axial function:
13 -> 12 => (0, -1, 1, 0)
14 -> 12 => (0, -1, 0, 1)
14 -> 13 => (0, 0, -1, 1)
23 -> 12 => (-1, 0, 1, 0)
23 -> 13 => (-1, 1, 0, 0)
24 -> 12 => (-1, 0, 0, 1)
24 -> 14 => (-1, 1, 0, 0)
24 -> 23 => (0, 0, -1, 1)
34 -> 13 => (-1, 0, 0, 1)
34 -> 14 => (-1, 0, 1, 0)
34 -> 23 => (0, -1, 0, 1)
34 -> 24 => (0, -1, 1, 0)

```
"""
function grassmannian(::Type{GKM_graph}, k::Int, n::Int)
  @req (k >= 0 && n >= k) "Dimension must be non-negative"
  
  return flag_variety(GKM_graph, [k, n-k])
end

@doc raw"""
    projective_space(::Type{gkm_graph}, d::Int)

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
```
"""
function projective_space(::Type{GKM_graph}, d::Int)
  @req d >= 0 "Dimension must be non-negative"
  
  return grassmannian(GKM_graph, 1, d+1)
end


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

@doc raw"""
    gkm_2d(w::Matrix{Int64}) -> AbstractGKM_graph

Return the 2-valent GKM cyclic connected GKM graph whose vertices are $1,2,\dots,n$ and whose edges
are $(1, 2), (2, 3), ..., (n, 1)$.
The weights of those edges are given by the rows of the matrix `w`.

# Example

The following example is the GKM graph from [GKZ22; Example 2.44, left figure](@cite), which cannot come from a Hamiltonian action.
One way of seeing this is that the combinatorial Betti numbers are not the geometric Betti numbers of any connected space.

```jldoctest
julia> G = gkm_2d([1 0; 0 1; -1 0; 0 -1; 1 0; 0 1; -1 0; 0 -1;])
GKM graph with 8 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (0, -1)
4 -> 3 => (1, 0)
5 -> 4 => (0, 1)
6 -> 5 => (-1, 0)
7 -> 6 => (0, -1)
8 -> 1 => (0, -1)
8 -> 7 => (1, 0)

julia> betti_numbers(G)
3-element Vector{Int64}:
 2
 4
 2
```
"""
function gkm_2d(w::Matrix{Int64})::AbstractGKM_graph
  s = size(w)
  n = s[1]
  @req n >= 3 "Need at least three vertices for 2d GKM graph"
  r = s[2]
  labels = ["$i" for i in 1:n]
  G = empty_gkm_graph(n, r, labels)

  for i in 1:n
    add_edge!(G, i, i % n + 1, G.M(w[i,:]))
  end
  return G
end

@doc raw"""
    gkm_3d_positive_non_toric(i::Int64) -> AbstractGKM_graph

Return the `i`-th GKM graph from [CK23; Appendix A](@cite)'s classification of 3-valent Hamiltonian positive GKM-graphs
with 2-dimensional torus-action that are not projections of GKM graphs coming from smooth projective polytopes.

The argument `i` runs from 1 ot 7 and is the index in the list.

# Example
We reproduce here the Betti numbers and the integrals $\int_M (c_1(M))^3$ as listed in [CK23; Appendix A](@cite).
```jldoctest
julia> for i in 1:7
           G = gkm_3d_positive_non_toric(i)
           integral = integrate(first_chern_class(G)^3, G)
           println("Graph $i: Betti numbers = $(betti_numbers(G)), [M] . (c_1(M))^3 = $integral")
         end
Graph 1: Betti numbers = [1, 1, 1, 1], [M] . (c_1(M))^3 = 54
Graph 2: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 30
Graph 3: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 40
Graph 4: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 46
Graph 5: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 48
Graph 6: Betti numbers = [1, 3, 3, 1], [M] . (c_1(M))^3 = 26
Graph 7: Betti numbers = [1, 3, 3, 1], [M] . (c_1(M))^3 = 38
```
"""
function gkm_3d_positive_non_toric(i::Int64)::AbstractGKM_graph
  if i == 1
    G = gkm_2d([1 1; -1 1; -1 -1; 1 -1])
    add_edge!(G, 1, 3, G.M([0, 1]))
    add_edge!(G, 4, 2, G.M([1, 0]))
    return G
  elseif i == 2
    G = gkm_2d([2 0; 0 1; -2 2; -1 0; 0 -2; 1 -1])
    add_edge!(G, 1, 4, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 3
    G = gkm_2d([2 0; 1 1; -1 1; -2 0; -1 -1; 1 -1])
    add_edge!(G, 1, 5, G.M([0, 1]))
    add_edge!(G, 2, 4, G.M([0, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 4
    G = gkm_2d([1 1; 0 1; -1 1; -1 0; -1 -1; 1 -1])
    add_edge!(G, 1, 4, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 5
    G = gkm_2d([1 0; 0 1; -1 1; -1 0; 0 -1; 1 -1])
    add_edge!(G, 1, 4, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 6
    G = gkm_2d([2 0; 1 1; 0 2; -1 1; -2 0; -1 -1; 0 -2; 1 -1])
    add_edge!(G, 1, 6, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([0, 1]))
    add_edge!(G, 3, 8, G.M([-1, 0]))
    add_edge!(G, 4, 7, G.M([-1, 0]))
    return G
  elseif i == 7
    G = gkm_2d([1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1])
    add_edge!(G, 1, 5, G.M([0, 1]))
    add_edge!(G, 2, 7, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 1]))
    add_edge!(G, 4, 8, G.M([-1, 0]))
    return G
  else
    @req false "Index must be between 1 and 7."
  end
end

@doc raw"""
    gkm_3d_twisted_flag() -> AbstractGKM_graph

Return the 3-valent GKM graph of the twisted flag varieties of Eschenburg, Tolman, and Woodward
(see [GKZ20; Example 4.8](@cite) and references therein).

# Example
Note that the resulting GKM graph does not occur in the output of `gkm_3d_positive_non_toric()` since one edge has
non-positive Chern number.

```jldoctest
julia> G = gkm_3d_twisted_flag()
GKM graph with 6 nodes, valency 3 and axial function:
2 -> 1 => (0, -1)
3 -> 2 => (1, 0)
4 -> 1 => (1, -2)
4 -> 3 => (-1, 1)
5 -> 2 => (1, -1)
5 -> 4 => (0, -1)
6 -> 1 => (1, -1)
6 -> 3 => (2, -1)
6 -> 5 => (1, 0)

julia> print_curve_classes(G)
2 -> 1: (0, 1), Chern number: 4
3 -> 2: (-1, 1), Chern number: 2
4 -> 1: (1, 0), Chern number: 2
4 -> 3: (-2, 1), Chern number: 0
5 -> 2: (1, 0), Chern number: 2
5 -> 4: (-1, 1), Chern number: 2
6 -> 1: (1, 1), Chern number: 6
6 -> 3: (1, 0), Chern number: 2
6 -> 5: (0, 1), Chern number: 4
```
"""
function gkm_3d_twisted_flag()::AbstractGKM_graph
  G = gkm_2d([0 1; -1 0; 1 -1; 0 1; -1 0; 1 -1])
  add_edge!(G, 1, 4, G.M([-1, 2]))
  add_edge!(G, 2, 5, G.M([-1, 1]))
  add_edge!(G, 3, 6, G.M([-2, 1]))
  return G
end

#TODO: implement [GKZ20, Prop. 4.5]'s classification of signed 3D GKM fibrations. 