function line_bundle(
  G::AbstractGKM_graph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}
)::GKM_vector_bundle where R<:GKM_weight_type
  return vector_bundle(G, M, GMtoM, reshape(weights, length(weights), 1))
end


@doc raw"""
    vector_bundle(G::AbstractGKM_graph, M::AbstractAlgebra.Generic.FreeModule{R}, GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}, weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}}; calculateConnection::Bool=true) -> GKM_vector_bundle

It constructs the vector bundle given by the following datam:#
# Arguments
- `g::G::AbstractGKM_graph`: A GKM graph
- `M::AbstractAlgebra.Generic.FreeModule{R}`: 

# Examples
```jldoctest
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> GMtoM = ModuleHomomorphism(G.M, M, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]])
GKM vector bundle of rank 1 over GKM graph with 3 nodes and valency 2 with weights:
1: (1, 0, 0, 0)
2: (0, 1, 0, 0)
3: (0, 0, 1, 0)

```
"""
function vector_bundle(
  G::AbstractGKM_graph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}};
  calculateConnection::Bool=true
)::GKM_vector_bundle where R<:GKM_weight_type

  @req domain(GMtoM) == G.M "GMtoM must go from gkm.M into M"
  @req rank(kernel(GMtoM)[1]) == 0 "GMtoM must be injective"
  @req codomain(GMtoM) == M "GMtoM must go from gkm.M into M"

  s = size(weights)
  nv = n_vertices(G.g)
  @req s[1] == nv "Weight matrix has wrong dimensions."
  
  for w in weights
    @req parent(w) == M "Weights need to live in M."
  end

  res =  GKM_vector_bundle(G, M, GMtoM, weights, nothing)
  # build connection if it is unique.
  if calculateConnection
    get_vector_bundle_connection(res)
  end
  return res
end

function vector_bundle_rank(V::GKM_vector_bundle)::Int64
  return size(V.w)[2]
end

function get_vector_bundle_connection(V::GKM_vector_bundle)
  if isnothing(V.con)
    V.con = _build_vector_bundle_connection(V)
  end
  return V.con
end

# Return the unique GKM conncetion of the vector bundle or nothing if it is not uniquely determined.
function _build_vector_bundle_connection(V::GKM_vector_bundle)

  con = Dict{Tuple{Edge, Int64}, Int64}()
  weights = V.w

  G = V.gkm
  rk = vector_bundle_rank(V)

  for e in edges(G.g)
    @req !is_zero(G.w[e]) "Weight zero edge found."
    v = src(e)
    w = dst(e)
    we = V.GMtoM(G.w[e])
    for i in 1:rk
      wi = weights[v, i]
      haveFoundJ = false
      for j in 1:rk
        wj = weights[w, j]
        wdif = wi - wj
        if rank(matrix([ wdif; we ])) == 1 # if true, (v,i) belongs to (w,j)
          if haveFoundJ
            # connection is not unique, so return nothing.
            return nothing
          else
            # have found a unique (so far) candidate for j.
            con[(e, i)] = j
            con[(reverse(e), j)] = i
            haveFoundJ = true
          end
        end
      end
      if !haveFoundJ
        return nothing
      end
    end
  end
  return con
end

function direct_sum(V::GKM_vector_bundle{R}...)::GKM_vector_bundle where R<:GKM_weight_type
  n = length(V)
  @req n >= 1 "Need at least one direct summand."
  for i in 1:n, j in 1:n
    @req V[i].gkm == V[j].gkm "Vector bundles need to have the same GKM base."
    @req V[i].M == V[j].M "Vector bundles need to have the same character lattice."
    @req V[i].GMtoM == V[j].GMtoM "V.GMtoM needs to be constant among direct summands."
  end
  G = V[1].gkm
  M = V[1].M
  GMtoM = V[1].GMtoM
  weights = hcat((V[i].w for i in 1:n)...)
  println(typeof(weights))
  res = vector_bundle(G, M, GMtoM, weights)

  # infer connection from direct summands
  if isnothing(get_vector_bundle_connection(res)) && !any([isnothing(get_vector_bundle_connection(V[i])) for i in 1:n])
    con = Dict{Tuple{Edge, Int64}, Int64}()
    offset = 0
    for i in 1:n
      conI = get_vector_bundle_connection(V[i])
      for k in keys(conI)
        e = k[1]
        a = k[2]
        b = conI[k]
        con[(e, a+offset)] = b+offset
        con[(reverse(e), b+offset)] = a+offset
      end
      offset += vector_bundle_rank(V[i])
    end
    res.con = con
  end
  return res
end


function Base.show(io::IO, V::GKM_vector_bundle)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM vector bundle")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM vector bundle of rank $(vector_bundle_rank(V)) over GKM graph with $(n_vertices(V.gkm.g)) vertices")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", V::GKM_vector_bundle)

  print(io, "GKM vector bundle of rank $(vector_bundle_rank(V)) over $(V.gkm) with weights:")
  rk = vector_bundle_rank(V)
  for v in 1:n_vertices(V.gkm.g)
    print(io, "\n$(V.gkm.labels[v]): ")
    for i in 1:rk
      print(io, V.w[v,i])
      if i<rk
        print(io, ", ")
      end
    end
  end
end

function dual_bundle(V::GKM_vector_bundle)::GKM_vector_bundle
  res = vector_bundle(V.gkm, V.M, V.GMtoM, -V.w; calculateConnection=false)
  res.con = V.con
  return res
end

function projective_bundle(V::GKM_vector_bundle)::AbstractGKM_graph
  con = get_vector_bundle_connection(V)
  @req !isnothing(con) "GKM vector bundle needs connection for projectivization."
  G = V.gkm
  nv = n_vertices(G.g)
  rk = vector_bundle_rank(V)

  Gres = Graph{Undirected}(nv * rk)
  labels = String[]
  sizehint!(labels, nv * rk)
  #build labels
  for v in 1:nv
    for i in 1:rk
      push!(labels, "[" * G.labels[v] * "]_$i")
    end
  end
  weightType = typeof(_get_weight_type(G))
  res = gkm_graph(Gres, labels, V.M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{weightType}}(); checkLabels=false)

  Gcon = get_GKM_connection(G)
  resCon = Dict{Tuple{Edge, Edge}, Edge}()

  # add edges corresponding to original edges
  for e in edges(G.g)
    v = src(e)
    w = dst(e)
    for i in 1:rk
      vInd = (v-1)*rk + i
      j = con[(e, i)]
      wInd = (w-1)*rk + j
      GKMadd_edge!(res, vInd, wInd, V.GMtoM(G.w[e]))

      if !isnothing(Gcon)
        eNew = Edge(vInd, wInd)
        for u in all_neighbors(G.g, v)
          uInd = (u-1)*rk + con[(Edge(v, u), i)]
          ei = Edge(vInd, uInd)
          up = dst(Gcon.con[(e, Edge(v, u))])
          upInd = (up - 1)*rk + con[(Edge(w, up), j)]
          epi = Edge(wInd, upInd)
          resCon[(eNew, ei)] = epi
          resCon[(reverse(eNew), epi)] = ei
        end
        for k in 1:rk
          i == k && continue
          kInd = (v-1)*rk + k
          l = con[(e, k)]
          lInd = (w-1)*rk + l
          ei = Edge(vInd, kInd)
          epi = Edge(wInd, lInd)
          resCon[(eNew, ei)] = epi
          resCon[(reverse(eNew), epi)] = ei
        end
      end
    end
  end
  # add edges over each vertex
  for v in 1:nv
    for i in 1:rk
      for j in (i+1):rk
        viInd = (v-1)*rk + i
        vjInd = (v-1)*rk + j
        wNew = V.w[v, j] - V.w[v, i]
        @req !iszero(wNew) "Vector bundle has two identical weights over vertex $v (indices $i, $j)"
        GKMadd_edge!(res, viInd, vjInd, wNew)

        if !isnothing(Gcon)
          e = Edge(viInd, vjInd)
          resCon[(e, e)] = reverse(e)
          resCon[(reverse(e), reverse(e))] = e
          for k in 1:rk
            (i == k || j == k) && continue
            vkInd = (v-1)*rk + k
            ei = Edge(viInd, vkInd)
            epi = Edge(vjInd, vkInd)
            resCon[(e, ei)] = epi
            resCon[(reverse(e), epi)] = ei
          end
          for w in all_neighbors(G.g, v)
            eDown = Edge(v, w)
            k = con[(eDown, i)]
            l = con[(eDown, j)]
            wkInd = (w-1)*rk + k
            wlInd = (w-1)*rk + l
            ei = Edge(viInd, wkInd)
            epi = Edge(vjInd, wlInd)
            resCon[(e, ei)] = epi
            resCon[(reverse(e), epi)] = ei
          end
        end
      end
    end
  end

  if !isnothing(Gcon)
    resConObj = build_GKM_connection(res, resCon)
    set_GKM_connection!(res, resConObj)
  end

  if !GKM_isValid(res)
    println("Warning: resulting projective bundle is not a valid GKM graph (see reason above).")
  end

  return res
end