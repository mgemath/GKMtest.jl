function *(G1::AbstractGKM_graph, G2::AbstractGKM_graph)
  
    n1 = n_vertices(G1.g)
    nv = n1 * n_vertices(G2.g)
  
    g = Graph{Undirected}(nv)
    M = free_module(ZZ, rank(G1.M)+rank(G2.M)) # direct_sum(G1.M, G2.M)
    f1 = hom(G1.M, M, [gens(M)[i] for i in 1:rank(G1.M)])
    f2 = hom(G2.M, M, [gens(M)[i + rank(G1.M)] for i in 1:rank(G2.M)])
    W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
    labels = Vector{String}(undef, nv)
  
    for e in edges(G1.g)
      v, w = src(e), dst(e)
      for _e in edges(G2.g)
        _v, _w = src(_e), dst(_e)
        
        V1 = v + (_v-1)*n1
        V2 = v + (_w-1)*n1
        V3 = w + (_v-1)*n1
        V4 = w + (_w-1)*n1
  
        add_edge!(g, V1, V2)
        add_edge!(g, V1, V3)
        
        add_edge!(g, V2, V4)
        add_edge!(g, V3, V4)
        
  
        W[Edge(max(V1, V3), min(V1, V3))] = f1(G1.w[e])
        W[Edge(max(V2, V4), min(V2, V4))] = f1(G1.w[e])
  
        W[Edge(max(V1, V2), min(V1, V2))] = f2(G2.w[_e])
        W[Edge(max(V3, V4), min(V3, V4))] = f2(G2.w[_e])
  
        labels[V1] = G1.labels[v]*","*G2.labels[_v]
        labels[V2] = G1.labels[v]*","*G2.labels[_w]
        labels[V3] = G1.labels[w]*","*G2.labels[_v]
        labels[V4] = G1.labels[w]*","*G2.labels[_w]
      end
    end
  
    return gkm_graph(g, labels, M, W)
  end