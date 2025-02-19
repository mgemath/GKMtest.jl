# # * It should probably take the GKM connection as input
# function integrateGKM(G::AbstractGKM_graph, classes::Vector{FreeModElem{QQMPolyRingElem}})
#   con = build_GKM_connection(G)
#   R = equivariant_cohomology_ring(G)

#   ########
#   # this part is needed for the generation of colorings
#   nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
#   for v in 1:n_vertices(G.g)
#     nc[v] = sort!(copy(all_neighbors(G.g, v)))
#   end
#   #########
#   DanielResult = 1 

#   res = zero(R.coeffRing)
#   println(res)

#   max_n_vert::Int64 = 3 #TODO find the maximal number of vertices

#   n_marks = length(classes)
#   for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
#     tree = LStoGraph(ls) # from level sequence to graph

#     CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
#     for col in CI   # colorings Iterator
#       vDict = Dict{Int, Int}([i for i in 1:n_vertices(tree)] .=> col) # TODO use directly col[i] instead of vDict[i]

#       for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
#         euler = zero(R.coeffRing)
#         for m in Combinatorics.multiset_permutations(m_inv, n_marks) #generation of all marks

#           tree_iso = count_iso(ls, col, m)  # isomorphism of the tree with coloring and marks

#           for edgeMult in all_ones_multip(G, tree)   # multiplicities

#             PROD = prod(e -> edgeMult[e], keys(edgeMult))
#             dt = decoratedTree(G, tree, vDict, edgeMult, m)

#             if euler == zero(R.coeffRing)
#               euler = Euler_inv(dt, R, con)//(PROD * tree_iso)
#             end
#             # res += euler
#             # println(euler)
            
#             # if col == (1, 4, 3) && m == [1, 1, 2]
#             #   DanielResult = euler
#             Class = R.coeffRing(1)
#               for j in 1:n_marks
#                 Class *= _ev(dt, j, classes[j])
#                 # println(Class)
#                 # res += Class*euler
#               end
#               res += Class*euler
#             # else
#             #   for j in 1:n_marks
#             #     Class = _ev(dt, j, classes[j])
#             #     println(Class)
#             #   end
#             # end 
#           end
#         end
#       end
#     end
#   end

#   println(res)
# end

# function integrateGKM(G::AbstractGKM_graph, n_marks::Int64, P_input, beta)
#   con = build_GKM_connection(G)
#   R = equivariant_cohomology_ring(G)
#   P = P_input.func

#   ########
#   # this part is needed for the generation of colorings
#   nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
#   for v in 1:n_vertices(G.g)
#     nc[v] = sort(all_neighbors(G.g, v))
#   end
#   #########

#   #########
#   # Dict in order to store H
#   h_dict::Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}} = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
#   ########
#   res = zero(R.coeffRing)
#   # println(res)

#   max_n_vert::Int64 = beta + 1 #TODO find the maximal number of vertices

#   # n_marks = length(classes)
#   for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
#     tree = LStoGraph(ls) # from level sequence to graph

#     CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
#     for col in CI   # colorings Iterator
#       top_aut::Int64 = count_iso(ls, col)
#       vDict = Dict{Int, Int}([i for i in 1:n_vertices(tree)] .=> col) # TODO use directly col[i] instead of vDict[i]

#       Multi = multi(G, tree, beta)
      

#       for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        
#         aut = count_iso(ls, col, m_inv)

#         for edgeMult_array in Multi # all_ones_multip(G, tree)   # multiplicities
#           PROD = prod(edgeMult_array)
#           euler = zero(R.coeffRing)

#           edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          

#           for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

            
#             dt = decoratedTree(G, tree, vDict, edgeMult, m, R)

#             Class = Base.invokelatest(P, dt)

#             Class == zero(R.coeffRing) && continue

#             if euler == zero(R.coeffRing)
#               euler = Euler_inv(dt, R, con)//(PROD * aut)
#               for e in edges(tree)
#                 triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
#                 if !haskey(h_dict, triple)
#                     h_dict[triple] = h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R, check=false)
#                 end
#                 euler *= h_dict[triple]
#               end
#             end
#             # println("ls = $(ls), col = $(col), aut = $(aut), PRODW = $(PROD), m=$(m), E = $(euler)")
#             # println(euler)
            
#             res += Class*euler
#           end
#         end
#       end
#     end
#   end
#   return res
# end

function integrateGKM(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input; show_bar::Bool = true)

  @req beta != zero(parent(beta)) "zero classes not yet supported. (TODO.)"

  H2 = GKM_second_homology(G)
  if !isEffectiveCurveClass(H2, beta)
    return 0
  end

  R = G.equivariantCohomology
  P = P_input.func
  con = get_GKM_connection(G)
  @req !isnothing(con) "GKM graph needs a connection!"

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########

  #########
  # Dict in order to store H
  h_dict::Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}} = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
  ########
  res = zero(R.coeffRing)
  # println(res)

  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  if show_bar #set up progress data
    number_trees = A000055(max_n_vert)
    # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
    # G.g, also counting ways to distribute the marked points among the vertices.
    # This does not cound the edge multiplicities.
    threshold = sum(vert -> number_trees[vert] * n_vertices(G.g) * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
end

  # n_marks = length(classes)
  # iterate undecorated trees:
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph
    tree_aut = count_iso(ls)

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)

      Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)

      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        
        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(R.coeffRing)

          edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          
          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

            
            dt = decoratedTree(G, tree, col, edgeMult, m)
            
            Class = Base.invokelatest(P, dt)

            Class == zero(R.coeffRing) && continue

            if euler == zero(R.coeffRing)
              euler = Euler_inv(dt)//(PROD * aut)
              for e in edges(tree)
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                if !haskey(h_dict, triple)
                    h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R, check=false)
                end
                euler *= h_dict[triple]
              end
            end
            # println("ls = $(ls), col = $(col), aut = $(aut), PRODW = $(PROD), m=$(m), E = $(euler)")
            # println(euler)
            
            res += Class*euler
          end
        end

        if show_bar #update the progress bar
          current_graph += tree_aut ÷ top_aut
          update!(progress_bar, current_graph,
              showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph)])
        end
      end
    end
  end
  return res
end