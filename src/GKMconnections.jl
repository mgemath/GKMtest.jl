"""
  Return the GKM_connection of the given GKM graph if it is unique.
  If it is not unique, return nothing.
"""
function build_GKM_connection(gkm::AbstractGKM_graph) :: GKM_connection
  
  if valency(gkm) >= 3
    @req is3_indep(gkm) "GKM graph has valency >= 3 is not 3-independent"
  elseif valency(gkm) == 2
    @req is2_indep(gkm) "GKM graph has valency 2 and is not 2-independent"
  end

  # assign to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  con = Dict{Edge, Dict{Edge, Edge}}()
  # w[e'_i] = w [e_i] - a_i * w[e]
  a = Dict{Edge, Dict{Edge, ZZRingElem}}()

  # iterate over all unoriented edges
  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eCon = Dict{Edge, Edge}()
    eA = Dict{Edge, ZZRingElem}()
    eRevCon = Dict{Edge, Edge}() # calculate connection for reversed edges too
    eRevA = Dict{Edge, ZZRingElem}()
    eW = gkm.w[e]

    # get all edges at src(e)
    for ei in [Edge(s1,v) for v in all_neighbors(gkm.g, s1)]

      # get all edges at dst(e), where epi stands for "e prime i"
      for epi in [Edge(s2, w) for w in all_neighbors(gkm.g, s2)]

        wdif = gkm.w[ei] - gkm.w[epi] # this will be a_i * w[e]
        
        if rank(matrix([ wdif; eW ])) == 1 # if true, epi belongs to ei.

          ai::ZZRingElem = ZZ(0) # divide wdiv by eW. Is this already implemented somewhere for FreeModuleElem{ZZRingElem}?
          for j in 1:rank(gkm.M)
            if eW[j] != 0
              ai = div(wdif[j],eW[j]) # type is again ZZRingElem
              break
            end
          end
          eCon[ei] = epi
          eRevCon[epi] = ei
          eA[ei] = ai
          eRevA[epi] = ai

          break
        end
      end
    end

    con[e] = eCon
    a[e] = eA
    con[reverse(e)] = eRevCon
    a[reverse(e)] = eRevA
  end

  return GKM_connection(gkm, con, a)
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