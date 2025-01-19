"""
Return the array betti such that betti[i+1] is the 2i-th (combinatorial) Betti number for i in {0,1,...,valency(gkm)}, as defined in [Guillemin--Zara, section 1.3].
Note:
 * i ranges from 0 to valency(gkm).
 * from [Guillemin--Zara, Theorem 1.3.2], the combinatorial Betti numbers equal the Betti numbers of the underlying
   GKM space if the torus action is Hamiltonian.
Warning:
 * betti[1] is the 0-th Betti number, since Julia arrays are 1-based and not 0-based.
""" 
function bettiNumbers(gkm::AbstractGKM_graph)

  for counter in 1:10^8 # arbitrary maximum number of attempts to avoid "while true"

    xi = ZZ.(rand(Int, n_vertices(gkm.g)))  #TODO: find something without using random numbers
    wxi = Dict{Edge, ZZRingElem}() # wxi stands for weight[e](xi)
    isPolarizing = true

    # calculate weight[e](xi) for all edges e
    for e in edges(gkm.g)

      wxi[e] = wxi[reverse(e)] = 0

      for j in 1:rank_torus(gkm)
        wxi[e] += xi[j] * gkm.w[e][j]
        wxi[reverse(e)] += xi[j] * gkm.w[reverse(e)][j]
      end

      if wxi[e] == 0 || wxi[reverse(e)] == 0
        isPolarizing = false
        break
      end
    end
    (!isPolarizing) && continue

    # from here on, xi is known to be polarizing.
    #betti = Dict{Int, Int}([i => 0 for i in 0:valency(gkm)])
    betti = zeros(Int, valency(gkm)+1)

    for v in vertices(gkm.g)

      i = count(w -> wxi[Edge(v,w)] < 0, all_neighbors(gkm.g, v) )
      betti[i+1] += 1
    end
    return betti
  end

end