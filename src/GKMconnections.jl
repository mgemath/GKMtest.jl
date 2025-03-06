"""
Return the GKM_connection of the given GKM graph if it is unique or has been set manually.
If it is unique and hasn't been calculated, it is saved in the gkm object.
If the connection is not unique and hasn't been defined manually, return nothing.
"""
function get_GKM_connection(gkm::AbstractGKM_graph)::Union{Nothing, GKM_connection}
  if isnothing(gkm.connection)
    if (valency(gkm) >= 3 && is3_indep(gkm)) || (valency(gkm) == 2 && is2_indep(gkm)) || (valency(gkm)==1)
      gkm.connection = _build_GKM_connection(gkm)
    end
  end
  return gkm.connection
end

function set_GKM_connection!(gkm::AbstractGKM_graph, con::GKM_connection)
  @req gkm == con.gkm "Connection belongs to the wrong GKM graph!"
  @req GKM_isValidConnection(con) "GKM connection is invalid!"
  
  gkm.connection = con
end

"""
Return the freshly calculated GKM_connection of the given GKM graph if it is unique.

Warning:
  1. This does not save the newly calculated GKM connection in the gkm object.
  2. If the connection is unique or was set before, one should instead use get_GKM_connection().
"""
function _build_GKM_connection(gkm::AbstractGKM_graph) :: GKM_connection
  
  if valency(gkm) >= 3
    @req is3_indep(gkm) "GKM graph has valency >= 3 is not 3-independent"
  elseif valency(gkm) == 2
    @req is2_indep(gkm) "GKM graph has valency 2 and is not 2-independent"
  end

  # assign to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  con = Dict{Tuple{Edge, Edge}, Edge}()

  # iterate over all unoriented edges
  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    # get all edges at src(e)
    for ei in [Edge(s1,v) for v in all_neighbors(gkm.g, s1)]

      # get all edges at dst(e), where epi stands for "e prime i"
      for epi in [Edge(s2, w) for w in all_neighbors(gkm.g, s2)]

        wdif = gkm.w[ei] - gkm.w[epi] # this will be a_i * w[e]
        
        if rank(matrix([ wdif; eW ])) == 1 # if true, epi belongs to ei.

          con[(e, ei)] = epi
          con[(reverse(e), epi)] = ei
          break
        end
      end
    end
  end

  return build_GKM_connection(gkm, con)
end

"""
Return the GKM connection object (including information of the ai's) determined by the given connection map.
Warning:
  1. This function does not check whether the given connection map is valid.
  2. This does not save the new connection to the gkm object.
"""
function build_GKM_connection(gkm::AbstractGKM_graph, con::Dict{Tuple{Edge, Edge}, Edge}) :: GKM_connection
  a = connection_a_from_con(gkm, con)
  return GKM_connection(gkm, con, a)
end

"""
Warning: this does not save the new connection to the gkm object.
"""
function build_GKM_connection(gkm::AbstractGKM_graph, a::Dict{Tuple{Edge, Edge}, ZZRingElem}) :: GKM_connection
  con = connection_map_from_a(gkm, a)
  return GKM_connection(gkm, con, a)
end

"""
Return the ai's belonging to the given GKM connection.
Warning: This function does not check whether the given connection map is valid.
"""
function connection_a_from_con(gkm::AbstractGKM_graph, con::Dict{Tuple{Edge, Edge}, Edge}; check::Bool = true)::Dict{Tuple{Edge, Edge}, ZZRingElem}

  a = Dict{Tuple{Edge, Edge}, ZZRingElem}()
  for e in edges(gkm.g)
    
    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    for ei in [Edge(s1,v) for v in all_neighbors(gkm.g, s1)]
      
      epi = con[(e, ei)]
      wdif = gkm.w[ei] - gkm.w[epi]

      if check
        @req rank(matrix([ wdif; eW ])) == 1 "connection is incompatible with GKM graph"
      end

      ai::ZZRingElem = ZZ(0)

      for j in 1:rank(gkm.M)
        if eW[j] != 0
          tmp = wdif[j] // eW[j]
          @req denominator(tmp) == 1 "GKM connection's a_i's must be integers!" # Assumption: x//y is integer if and only if denominator(x//y) == 1 in Oscar.
          ai = ZZ(tmp)
          break
        end
      end

      a[(e, ei)] = ai
      a[(reverse(e), epi)] = ai
    end
  end
  return a
end

"""
Build the connection map from the given collection of a's [cf. Liu--Sheshmani 2.(b) on p.4]
Warning: The returned value is only unique if the GKM has no repeated weights at any vertex (which is required for it to be valid).
"""
function connection_map_from_a(gkm::AbstractGKM_graph, a::Dict{Tuple{Edge, Edge}, ZZRingElem})::Dict{Tuple{Edge, Edge}, Edge}

  con = Dict{Tuple{Edge, Edge}, Edge}()
  for e in edges(gkm.g)
    for ei in [Edge(src(e),v) for v in all_neighbors(gkm.g, src(e))]

      ai = a[(e, ei)]
      wEpi = gkm.w[ei] - ai * gkm.w[e] # following [Liu--Sheshmani 2.(b) on p.4]

      resultFound = false

      for epi in [Edge(dst(e),v) for v in all_neighbors(gkm.g, dst(e))]

        if wEpi == gkm.w[epi]

          con[(e, ei)] = epi
          con[(reverse(e), epi)] = ei
          resultFound = true
          break
        end
      end
      @req resultFound "No edge found for ($e,$ei) using connection a's"
    end
  end
  return con
end

"""
Return true if the given GKM connection is valid. This holds if and only if all of the following hold:
  1. con and a are set for all (Edge(v,w), Edge(v,u)) where vw and vu are connected in the graph
  2. con maps every (e,e) to reverse(e)
  3. a maps every (e,e) to 2
  4. Every pair (e,ei) with same source satisfies the relation of the associated a's, i.e. w[ei'] = w[ei] - a[(e,ei)] * w[e]
"""
function GKM_isValidConnection(con::GKM_connection; printDiagnostics::Bool=true)::Bool

  for e in edges(con.gkm.g)
    if !haskey(con.con, (e,e))
      printDiagnostics && println("Connection misses key (e,e) for e=$e.")
      return false
    elseif !haskey(con.a, (e,e))
      printDiagnostics && println("Connection misses a(e,e) for e=$e.")
      return false
    elseif !haskey(con.con, (reverse(e),reverse(e)))
        printDiagnostics && println("Connection misses key (e,e) for e=$(reverse(e)).")
        return false
    elseif !haskey(con.a, (reverse(e),reverse(e)))
        printDiagnostics && println("Connection misses a(e,e) for e=$(reverse(e)).")
        return false
    elseif con.con[(e,e)] != reverse(e)
      printDiagnostics && println("Connection doesn't map (e,e) to reverse(e) for e=$e.")
      return false
    elseif con.con[(reverse(e),reverse(e))] != e
      printDiagnostics && println("Connection doesn't map (e,e) to reverse(e) for e=$e.")
      return false
    elseif con.a[(e,e)] != ZZ(2)
      printDiagnostics && println("Connection does not satisfy a(e,e)=2 for e=$e.")
      return false
    elseif con.a[(reverse(e), reverse(e))] != ZZ(2)
      printDiagnostics && println("Connection does not satisfy a(e,e)=2 for e=$(reverse(e)).")
      return false 
    end
  end

  for v in 1:n_vertices(con.gkm.g)
    for w in 1:n_vertices(con.gkm.g)
      (v == w) && continue
      e = Edge(v,w)
      if has_edge(con.gkm.g, e)
        for u in all_neighbors(con.gkm.g, v)

          ei = Edge(v, u)
          if !haskey(con.con, (e, ei))
            printDiagnostics && println("Connection map misses value for ($e, $ei).")
            return false
          elseif !haskey(con.a, (e, ei))
            printDiagnostics && println("Connection misses value for a($e, $ei).")
            return false
          end
          epi = con.con[(e,ei)]
          ai = con.a[(e,ei)]
          if con.gkm.w[epi] != con.gkm.w[ei] - base_ring(con.gkm.M)(ai) * con.gkm.w[e]
            printDiagnostics && println("Connection map and a(e,ei) is inconsistent for (e, ei)=($e, $ei).")
            return false
          end
        end
      end
    end
  end
  return true
end

function Base.show(io::IO, con::GKM_connection)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM connection")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM connection for GKM graph with $(n_vertices(con.gkm.g)) nodes and valency $(valency(con.gkm))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", con::GKM_connection)

  print(io, "GKM connection for GKM graph with $(n_vertices(con.gkm.g)) nodes and valency $(valency(con.gkm)):\n")
  print(io, "Connection:\n")
  show(io, MIME"text/plain"(), con.con)
  print(io, "\n a_i's:\n")
  show(io, MIME"text/plain"(), con.a)
  
end