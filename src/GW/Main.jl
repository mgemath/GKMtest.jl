function integrateGKM(G::AbstractGKM_graph, classes::Vector{FreeModElem{QQMPolyRingElem}})
  con = build_GKM_connection(G)
  R = equivariant_cohomology_ring(G)

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort!(copy(all_neighbors(G.g, v)))
  end
  #########
  DanielResult = 1 

  max_n_vert::Int64 = 3 #TODO find the maximal number of vertices

  n_marks = length(classes)
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    for col in CI   # colorings Iterator
      vDict = Dict{Int, Int}([i for i in 1:n_vertices(tree)] .=> col) # TODO use directly col[i] instead of vDict[i]

      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        for m in Combinatorics.multiset_permutations(m_inv, n_marks) #generation of all marks

          tree_iso = count_iso(ls, col, m)  # isomorphism of the tree with coloring and marks

          for edgeMult in all_ones_multip(G, tree)   # multiplicities

            PROD = prod(e -> edgeMult[e], keys(edgeMult))
            dt = decoratedTree(G, tree, vDict, edgeMult, m)

            euler = Euler_inv(dt, R, con)#//(PROD * tree_iso)
            # println(euler)
            
            if col == (1, 4, 3) && m == [1, 1, 2]
              DanielResult = euler
              for j in 1:n_marks
                Class = _ev(dt, j, classes[j])
                println(Class)
                DanielResult *= Class
              end
            else
              for j in 1:n_marks
                Class = _ev(dt, j, classes[j])
                println(Class)
              end
            end 
          end
        end
      end
    end
  end

  println("Daniel example is:\n$DanielResult")
end

function integrateGKM(G::AbstractGKM_graph, n_marks::Int64, P_input)
  con = build_GKM_connection(G)
  R = equivariant_cohomology_ring(G)
  P = P_input.func

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort!(copy(all_neighbors(G.g, v)))
  end
  #########
  DanielResult = 1 

  max_n_vert::Int64 = 3 #TODO find the maximal number of vertices

  # n_marks = length(classes)
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    for col in CI   # colorings Iterator
      vDict = Dict{Int, Int}([i for i in 1:n_vertices(tree)] .=> col) # TODO use directly col[i] instead of vDict[i]

      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        for m in Combinatorics.multiset_permutations(m_inv, n_marks) #generation of all marks

          tree_iso = count_iso(ls, col, m)  # isomorphism of the tree with coloring and marks

          for edgeMult in all_ones_multip(G, tree)   # multiplicities

            PROD = prod(e -> edgeMult[e], keys(edgeMult))
            dt = decoratedTree(G, tree, vDict, edgeMult, m)

            euler = Euler_inv(dt, R, con)#//(PROD * tree_iso)
            # println(euler)
            
            Class = Base.invokelatest(P, dt)
            # println(Class)
            if col == (1, 4, 3)
              if m == [1, 1, 2]
                println("euler: \n", euler)
                println("class: \n", Class)
                println("Product: \n", euler*Class)
                # DanielResult = euler
              end
            end  
          end
        end
      end
    end
  end

  # println("Daniel example is:\n$DanielResult")
end