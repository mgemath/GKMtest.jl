@doc raw"""
    betti_numbers(G::AbstractGKM_graph) -> Vector{Int64}

Return the array `betti` such that `betti_numbers[i+1]` is the 2i-th (combinatorial) Betti number for i from 0 to the valency of `G`, as defined in [GZ01; section 1.3](@cite).
!!! note
    * i ranges from 0 to the valency of `G`, that can be obtained by `valency(G)`.
    * from [GZ01; Theorem 1.3.2](@cite), the combinatorial Betti numbers equal the Betti numbers of the underlying GKM space if the torus action is Hamiltonian.

!!! warning
    `betti_numbers[1]` is the 0-th Betti number, since Julia arrays are 1-based and not 0-based.

""" 
function Oscar.betti_numbers(G::AbstractGKM_graph)::Vector{Int64}

  for counter in 1:10^8 # arbitrary maximum number of attempts to avoid "while true"

    xi = ZZ.(rand(Int, rank_torus(G)))  #TODO: find something without using random numbers
    wxi = Dict{Edge, ZZRingElem}() # wxi stands for weight[e](xi)
    isPolarizing = true

    # calculate weight[e](xi) for all edges e
    for e in edges(G.g)

      wxi[e] = wxi[reverse(e)] = 0

      for j in 1:rank_torus(G)
        wxi[e] += xi[j] * G.w[e][j]
        wxi[reverse(e)] += xi[j] * G.w[reverse(e)][j]
      end

      if wxi[e] == 0 || wxi[reverse(e)] == 0
        isPolarizing = false
        break
      end
    end
    (!isPolarizing) && continue

    # from here on, xi is known to be polarizing.
    #betti = Dict{Int, Int}([i => 0 for i in 0:valency(G)])
    betti = zeros(Int, valency(G)+1)

    for v in 1:n_vertices(G.g)

      i = count(w -> wxi[Edge(v,w)] < 0, all_neighbors(G.g, v) )
      betti[i+1] += 1
    end
    return betti
  end

end