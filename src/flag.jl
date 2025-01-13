# compute the flag variety whose quotient dimensions are expressed by the array s
function flag_gkm_graph(s::Vector{Int64})

  @req !isempty(s) "the vector of dimensions is empty"
  @req all(i -> s[i] > 0, eachindex(s)) "all dimensions must be positive"


    # accepted_perm = Vector{Int64}[]
  K = [sum(s[1:i]) for i in 0:length(s)]
  d = Dict{Int64, Vector{Int64}}()
  index = 1
  
  for c in Combinatorics.permutations(1:K[end])
    if all(i-> issorted(c[(K[i]+1):K[i+1]]), 1:(length(K)-1))
        # push!(accepted_perm, c)
      d[index] = c
      index += 1
    end
  end
      
    # nv = length(accepted_perm)
  nv = length(keys(d))
    # d = Dict(1:nv  .=>  accepted_perm)
  g = Graph{Undirected}(nv)
  
    # for v in 1:nv
    #   for w in (v+1):nv
    #     if all(i-> count(j -> (j in d[w][(K[i]+1):K[i+1]]), d[v][(K[i]+1):K[i+1]]) >= K[i+1]-(K[i]+1), 1:(length(K)-1)) # issorted(c[(K[i]+1):K[i+1]]), 1:(length(K)-1)) #count(i -> (i in d[w]), d[v]) == (k-1)
    #       if count(i-> (count(j -> (j in d[w][(K[1]+1):K[i+1]]), d[v][(K[1]+1):K[i+1]]) == K[i+1]-(K[1]+1)), 1:(length(K)-1)) < 2
    #         add_edge!(g, v, w)
    #       end
    #     end
    #   end
    # end

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

    # labels = [prod(i -> c[i]>9 ? "|$(c[i])|" : "$(c[i])", 1:K[end-1]) for c in accepted_perm]
  labels = [prod(i -> d[n][i]>9 ? "|$(d[n][i])|" : "$(d[n][i])", 1:K[end-1]) for n in 1:nv]

    # for e in edges(g)
    #     v = src(e); w = dst(e);
    #     for i in 1:(length(K)-1)
    #         if (count(j -> (j in d[w][(K[1]+1):K[i+1]]), d[v][(K[1]+1):K[i+1]]) == K[i+1]-(K[1]+1))
    #             value1 = setdiff(d[v][(K[1]+1):K[i+1]], d[w][(K[1]+1):K[i+1]])[1]
    #             value2 = setdiff(d[w][(K[1]+1):K[i+1]], d[v][(K[1]+1):K[i+1]])[1]
    #             W[e] = gens(M)[value1] - gens(M)[value2]
    #         end
    #     end
    # end
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
